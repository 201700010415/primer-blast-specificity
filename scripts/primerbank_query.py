# -*- coding: utf-8 -*-
"""按基因（GeneID 或基因名）在线查询 PrimerBank 引物对。

流程: 基因名 --(本地映射表 gene_map.tsv)--> GeneID --(PrimerBank 在线接口)--> 全部引物对。
引物结果不缓存，每次都实时查；映射表是本地的（刷新用 build_gene_map.py）。

用法:
  python primerbank_query.py --gene HMOX1 --species mouse
  python primerbank_query.py --gene 15368 --species mouse --csv out.csv
"""
import argparse
import csv
import gzip
import io
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

SEARCH_URL = 'https://pga.mgh.harvard.edu/cgi-bin/primerbank/new_search2.cgi'
UA = {'User-Agent': 'Mozilla/5.0 (primer-blast-specificity skill)'}
SPECIES = {'human': ('Human', '9606'), 'mouse': ('Mouse', '10090')}
GENE_MAP = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'gene_map.tsv')


def resolve_gene(query, tax_id):
    """本地映射表解析: GeneID 或基因名/别名 -> (gene_id, symbol, description)。
    找到多个候选时抛 ValueError 列出候选。"""
    sym2id, id2row = {}, {}
    with open(GENE_MAP, encoding='utf-8') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row['tax_id'] != tax_id:
                continue
            id2row[row['gene_id']] = row
            keys = [row['symbol']] + [s for s in row['synonyms'].split('|') if s and s != '-']
            for k in keys:
                sym2id.setdefault(k.upper(), set()).add(row['gene_id'])

    q = query.strip()
    if re.fullmatch(r'\d+', q):                      # 输入就是 GeneID
        if q not in id2row:
            raise ValueError('GeneID %s 不在 %s 物种映射表中（先重跑 build_gene_map.py）' % (q, tax_id))
        row = id2row[q]
        return q, row['symbol'], row['description']

    hits = sorted(sym2id.get(q.upper(), ()), key=int)
    if not hits:
        raise ValueError('基因名 "%s" 在映射表中无匹配（支持官方 symbol 和别名）；'
                         '确认拼写或直接给 GeneID' % q)
    if len(hits) > 1:
        cand = ', '.join('%s(%s)' % (g, id2row[g]['symbol']) for g in hits[:10])
        raise ValueError('基因名 "%s" 匹配到 %d 个 GeneID，请改用 --gene 指定 GeneID: %s%s'
                         % (q, len(hits), cand, ' ...' if len(hits) > 10 else ''))
    g = hits[0]
    return g, id2row[g]['symbol'], id2row[g]['description']


def http_post(url, data, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=UA)
            return urllib.request.urlopen(req, timeout=90).read().decode('utf-8', 'ignore')
        except Exception as e:
            last = e
            if i < tries - 1:
                print('[info] 网络异常，%ds 后重试: %s' % (3 * (i + 1), last), file=sys.stderr)
                time.sleep(3 * (i + 1))
    raise last


def search_primerbank(symbol, species_label):
    """按官方 symbol 提交（PrimerBank 的 GeneID 搜索方式服务端已失效）。"""
    payload = urllib.parse.urlencode({
        'selectBox': 'NCBI Gene Symbol', 'species': species_label,
        'searchBox': symbol, 'Submit': 'Submit'}).encode()
    return http_post(SEARCH_URL, payload)


def parse_result(html):
    """返回 (gene_meta: dict, pairs: list[dict])。"""
    txt = re.sub(r'<[^>]+>', ' ', html)
    txt = txt.replace('&nbsp;', ' ').replace('-&gt;', '->').replace('&gt;', '>')
    txt = re.sub(r'\s+', ' ', txt)

    meta = {}
    m = re.search(r'Gene Descriptions:.*?GeneID\s+(\d+).*?GenBank Accession\s+(\S+)', txt)
    if m:
        meta['gene_id'], meta['accession'] = m.group(1), m.group(2)
    m = re.search(r'Gene Description (.*?) Species ', txt)
    if m:
        meta['gene_desc'] = m.group(1).strip()

    if 'No matching' in txt or 'no match' in txt.lower():
        return meta, []

    pairs = []
    blocks = re.split(r'Primer Pair \d+\s*(?:\(Click here[^)]*\))?\s*:', txt)
    heads = re.findall(r'Primer Pair (\d+)', txt)
    for idx, block in enumerate(blocks[1:]):
        p = {}
        p['pair_no'] = heads[idx] if idx < len(heads) else str(idx + 1)
        m = re.search(r'PrimerBank ID (\S+)', block)
        p['pb_id'] = m.group(1) if m else '?'
        m = re.search(r'Amplicon Size (\d+)', block)
        p['amplicon'] = m.group(1) if m else '?'
        m = re.search(r'Forward Primer\s+([ACGTUacgtu]+)\s+(\d+)\s+([\d.]+)\s+([\d-]+)', block)
        if m:
            p['fwd'], p['fwd_len'], p['fwd_tm'], p['fwd_loc'] = m.group(1), m.group(2), m.group(3), m.group(4)
        m = re.search(r'Reverse Primer\s+([ACGTUacgtu]+)\s+(\d+)\s+([\d.]+)\s+([\d-]+)', block)
        if m:
            p['rev'], p['rev_len'], p['rev_tm'], p['rev_loc'] = m.group(1), m.group(2), m.group(3), m.group(4)
        p['validated'] = 'Validation Results' in block
        if p.get('fwd') and p.get('rev'):
            pairs.append(p)
        else:
            print('[warn] 第 %s 对解析不完整，已跳过' % p['pair_no'], file=sys.stderr)
    return meta, pairs


def main():
    p = argparse.ArgumentParser(description='PrimerBank 在线查询（GeneID 优先）')
    p.add_argument('--gene', required=True, help='基因名（官方 symbol 或别名）或纯数字 GeneID')
    p.add_argument('--species', required=True, choices=['human', 'mouse'])
    p.add_argument('--csv', default='', help='把引物对导出为 blast 用 CSV（fwd,rev,organism,target_gene）')
    p.add_argument('--blast', action='store_true',
                   help='查到引物后自动逐对跑 Primer-BLAST 特异性并出合并存档报告')
    p.add_argument('--delay', type=float, default=5, help='--blast 时每对提交间隔秒数')
    p.add_argument('--report', default='', help='--blast 时合并存档报告路径（默认按时间命名）')
    args = p.parse_args()

    species_label, tax_id = SPECIES[args.species]
    organism = 'Homo sapiens' if tax_id == '9606' else 'Mus musculus'

    # 映射表过期提醒（gene_info 基因集会随时间变动）
    age_days = (time.time() - os.path.getmtime(GENE_MAP)) / 86400.0
    if age_days > 180:
        print('[info] gene_map.tsv 已 %d 天未更新，新/改名基因可能查不到，'
              '建议重跑 build_gene_map.py' % age_days, file=sys.stderr)

    try:
        gene_id, symbol, desc = resolve_gene(args.gene, tax_id)
    except ValueError as e:
        print('基因解析失败: %s' % e, file=sys.stderr)
        sys.exit(2)
    # 注: PrimerBank 的 "NCBI Gene ID" 搜索方式服务端已失效（返回空页，2026-10 实测），
    # 因此用映射表解析出的官方 symbol 提交查询；结果页显示的 GeneID 可与本地解析互相核对。
    print('[info] %s -> GeneID %s (%s) [%s]，按官方 symbol 提交查询'
          % (args.gene, gene_id, symbol, organism), file=sys.stderr)

    html = search_primerbank(symbol, species_label)
    meta, pairs = parse_result(html)
    if meta.get('gene_id') and meta['gene_id'] != gene_id:
        print('[warn] PrimerBank 返回 GeneID %s != 本地解析 %s（symbol 可能对应多基因，'
              '请以结果页为准人工核对）' % (meta['gene_id'], gene_id), file=sys.stderr)

    lines = ['## PrimerBank 查询结果: %s (%s)\n' % (symbol, organism),
             '- 查询路径: `%s` --(本地映射表)--> GeneID %s --(官方 symbol 提交)--> PrimerBank'
             % (args.gene, gene_id),
             '- GenBank: %s | 描述: %s' % (
                 meta.get('accession', '?'), meta.get('gene_desc', desc)),
             '- 共 %d 对引物（来源 PrimerBank，实时查询）\n' % len(pairs)]
    if not pairs:
        lines.append('PrimerBank 未收录该基因的引物。')
    else:
        lines.append('| # | PrimerBank ID | 验证 | 产物bp | 正向引物 | 反向引物 | F Tm | R Tm | F/R 位置 |')
        lines.append('|---|---|---|---|---|---|---|---|---|')
        for p_ in pairs:
            lines.append('| %s | %s | %s | %s | `%s` | `%s` | %s | %s | %s/%s |' % (
                p_['pair_no'], p_['pb_id'], '是' if p_['validated'] else '',
                p_['amplicon'], p_['fwd'], p_['rev'],
                p_['fwd_tm'], p_['rev_tm'], p_['fwd_loc'], p_['rev_loc']))

    report = '\n'.join(lines) + '\n'
    out_md = '%s_%s.md' % (symbol, organism.replace(' ', '_'))
    with open(out_md, 'w', encoding='utf-8') as f:
        f.write(report)
    sys.stdout.write(report)
    print('\n报告已保存: %s' % os.path.abspath(out_md))

    if args.csv:
        with open(args.csv, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['fwd', 'rev', 'organism', 'target_gene'])
            for p_ in pairs:
                w.writerow([p_['fwd'], p_['rev'], organism, symbol])
        print('\n已导出 blast 用 CSV: %s（%d 行）' % (os.path.abspath(args.csv), len(pairs)))

    if args.blast:
        if not pairs:
            print('\n无引物可跑，跳过特异性验证', file=sys.stderr)
        else:
            print('\n[info] Primer-BLAST 特异性验证: %d 对（每对约 30-120s）'
                  % len(pairs), file=sys.stderr)
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import blast_primer as bp
            results = []
            for i, p_ in enumerate(pairs):
                print('\n--- [%d/%d] PrimerBank %s ---'
                      % (i + 1, len(pairs), p_['pb_id']))
                r = bp.run_one(dict(fwd=p_['fwd'], rev=p_['rev'],
                                    organism=organism, target_gene=symbol))
                results.append(r)
                bp.print_md(r)
                if i < len(pairs) - 1:
                    time.sleep(args.delay)
            ok = [r for r in results if not r.get('error')]
            if ok:
                report = args.report or os.path.abspath(
                    'PrimerBlast_report_%s_%s_%s.html' % (
                        symbol, organism.replace(' ', '_'),
                        time.strftime('%Y%m%d_%H%M%S')))
                bp.render_html(ok, report)
                print('\n合并存档报告: %s' % os.path.abspath(report))
            print('特异性验证完成: %d 成功 / %d 失败（共 %d）'
                  % (len(ok), len(results) - len(ok), len(results)))


if __name__ == '__main__':
    main()
