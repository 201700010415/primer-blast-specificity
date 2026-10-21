# -*- coding: utf-8 -*-
"""NCBI Primer-BLAST 引物特异性检查脚本。

输入正/反向引物与物种，提交 Primer-BLAST，轮询取回最终结果页，
解析出每个基因的命中产物与 F/R 引物错配数，输出摘要 + 原始结果页 URL。

仅依赖 Python 3 标准库。
用法示例（单对）：
  python blast_primer.py --fwd GGAGCGAGATCCCTCCAAAAT --rev GGCTGTTGTCATACTTCTCATGG \
      --organism "Homo sapiens" [--db refseq_mrna] [--target-gene GAPDH] \
      [--product-min 70] [--product-max 1000] [--mismatch 2] [--mismatch-3end 2] \
      [--out result.html]
批量模式（表头 fwd,rev,organism[,target_gene]，输出合并存档 HTML 报告）：
  python blast_primer.py --csv primers.csv [--report report.html] [--delay 5]
"""
import argparse
import csv
import http.client
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from html import escape as hesc

SUBMIT_URL = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'
FORM_URL = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/'
TAXID_URL = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
             '?db=taxonomy&term={term}&retmode=json')
UA = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) primer-blast-specificity-skill'

_ctx = None


def _ssl_ctx():
    global _ctx
    if _ctx is None:
        try:
            c = ssl.create_default_context()
            urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/', timeout=15, context=c).read(100)
            _ctx = c
        except Exception:
            c = ssl.create_default_context()
            c.check_hostname = False
            c.verify_mode = ssl.CERT_NONE
            _ctx = c
    return _ctx


def _retry(fn, label, tries=3, wait=5):
    """瞬时网络错误（NCBI 偶发超时/断连）自动重试。"""
    for i in range(tries):
        try:
            return fn()
        except urllib.error.HTTPError as e:
            if e.code in (429,) or e.code >= 500:
                last = e
            else:
                raise
        except (urllib.error.URLError, OSError, http.client.HTTPException) as e:
            last = e
        if i == tries - 1:
            raise last
        print('[info] %s网络异常，%ds 后重试(%d/%d): %s'
              % (label, wait * (i + 1), i + 1, tries, last), file=sys.stderr)
        time.sleep(wait * (i + 1))


def http_get(url, timeout=120):
    def _do():
        req = urllib.request.Request(url, headers={'User-Agent': UA, 'Referer': FORM_URL})
        r = urllib.request.urlopen(req, timeout=timeout, context=_ssl_ctx())
        return r.geturl(), r.read().decode('utf-8', 'ignore')
    return _retry(_do, 'GET ')


def http_post(url, data, timeout=420):
    def _do():
        req = urllib.request.Request(url, data=data, method='POST', headers={
            'Content-Type': 'multipart/form-data; boundary=' + _boundary,
            'User-Agent': UA, 'Referer': FORM_URL})
        r = urllib.request.urlopen(req, timeout=timeout, context=_ssl_ctx())
        return r.geturl(), r.read().decode('utf-8', 'ignore')
    return _retry(_do, 'POST ')


_boundary = ''


def resolve_organism(name):
    """把物种名解析成 'Name (taxid:N)' 格式；已带 taxid 或解析失败则原样返回。"""
    if 'taxid:' in name.lower():
        return name
    try:
        url = TAXID_URL.format(term=urllib.parse.quote(name))
        _, body = http_get(url, timeout=30)
        ids = re.search(r'"idlist"\s*:\s*\[\s*"(\d+)"', body)
        if ids:
            return '%s (taxid:%s)' % (name, ids.group(1))
    except Exception as e:
        print('[warn] taxid 解析失败，按原名提交: %s' % e, file=sys.stderr)
    return name


def build_form(job, organism):
    fields = {
        'CMD': 'request',
        'PRIMER_LEFT_INPUT': job['fwd'].strip().upper().replace(' ', ''),
        'PRIMER_RIGHT_INPUT': job['rev'].strip().upper().replace(' ', ''),
        'PRIMER_PRODUCT_MIN': str(job.get('product_min', 70)),
        'PRIMER_PRODUCT_MAX': str(job.get('product_max', 1000)),
        'PRIMER_NUM_RETURN': str(job.get('num_return', 10)),
        'SEARCH_SPECIFIC_PRIMER': 'on',
        'SEARCHMODE': '0',
        'PRIMER_SPECIFICITY_DATABASE': job.get('db', 'refseq_mrna'),
        'ORGANISM': organism,
        'ENTREZ_QUERY': job.get('entrez_query') or '',
        'TOTAL_PRIMER_SPECIFICITY_MISMATCH': str(job.get('mismatch', 2)),
        'PRIMER_3END_SPECIFICITY_MISMATCH': str(job.get('mismatch_3end', 2)),
        'MISMATCH_REGION_LENGTH': '5',
        'TOTAL_MISMATCH_IGNORE': '6',
        'MAX_TARGET_SIZE': '4000',
        'NUM_TARGETS': '20',
        'LINK_LOC': '',
        # NCBI 对自动化请求的礼仪标识（NCBI_ToU）：email 设环境变量 NCBI_EMAIL 后自动附带
        'tool': 'primer-blast-specificity-skill',
    }
    email = os.environ.get('NCBI_EMAIL', '').strip()
    if email:
        fields['email'] = email
    return fields


def submit(fields):
    global _boundary
    _boundary = '----pb' + uuid.uuid4().hex
    body = b''.join(
        ('--%s\r\nContent-Disposition: form-data; name="%s"\r\n\r\n%s\r\n'
         % (_boundary, k, v)).encode('utf-8')
        for k, v in fields.items()) + ('--%s--\r\n' % _boundary).encode()
    return http_post(SUBMIT_URL, body)


def find_refresh(html):
    """提取 meta refresh 目标 URL（任务轮询地址）；无则返回 None。"""
    m = re.search(r'<meta[^>]*http-equiv\s*=\s*Refresh[^>]*content\s*=\s*"\s*\d+\s*;\s*URL\s*=\s*([^"]+)"',
                  html, re.I)
    if not m:
        return None
    return urllib.parse.unquote(m.group(1).strip())


def wait_result(job_url, max_min=10):
    t0 = time.time()
    url = job_url
    retries = 0
    while time.time() - t0 < max_min * 60:
        try:
            final_url, page = http_get(url, timeout=120)
        except Exception as e:
            print('[info] 轮询异常，15s 后重试: %s' % e, file=sys.stderr)
            time.sleep(15)
            continue
        if 'Search Summary' in page and re.search(r'<title>\s*Primer-Blast results\s*</title>', page, re.I):
            return final_url, page
        nxt = find_refresh(page)
        if nxt and 'job_key' in nxt:
            url = nxt
            m = re.search(r'CONTENT="\s*(\d+)\s*;', page, re.I)
            wait_s = max(10, min(int(m.group(1)) if m else 20, 30))
            print('[info] 任务运行中，%ds 后轮询 (%.0fs)' % (wait_s, time.time() - t0), file=sys.stderr)
            time.sleep(wait_s)
            continue
        if nxt:  # 跳到非 job_key 页面（偶发），继续查原任务地址
            retries += 1
            if retries > 6:
                return final_url, page
            time.sleep(15)
            continue
        # 无 refresh 又不满足完成条件：可能异常页
        retries += 1
        if retries > 6:
            return final_url, page
        time.sleep(15)
    raise RuntimeError('轮询超时（%d 分钟）' % max_min)


def extract_gene(desc, acc):
    """从 'Homo sapiens actin beta (ACTB), mRNA' 提取基因符号。"""
    syms = re.findall(r'\(([A-Za-z0-9][A-Za-z0-9\-\.]{1,30})\)', desc)
    if syms:
        return syms[-1]
    return acc.split('.')[0]


def count_mismatch(alignment):
    """比对串中 '.' 为匹配，字母/('-') 计为错配。"""
    return len(re.findall(r'[A-Za-z\-]', alignment))


def parse_report(html):
    """解析结果页，返回 (摘要行, 分区命中列表)。"""
    summary_line = ''
    m = re.search(r'(No target templates were found[^<]*|Target templates were found[^<]*)', html)
    if m:
        summary_line = re.sub(r'\s+', ' ', m.group(1)).strip()

    hits = []
    for sec in re.finditer(
            r'<div class="prPairTl">([^<]+)</div>\s*<hr/>\s*<div class="prPairDtl">(.*?)</div>',
            html, re.S):
        section, inner = sec.group(1).strip(), sec.group(2)
        if not inner.strip():
            continue
        for h in re.finditer(
                r'><a[^>]*>([^<]+)</a>\s*([^<]*?)\s*<pre>(.*?)</pre>', inner, re.S):
            acc, desc, pre = h.group(1).strip(), h.group(2).strip(), h.group(3)
            plen = re.search(r'product length\s*=\s*(\d+)', pre)
            aln = {'fwd': [], 'rev': []}
            for a in re.finditer(
                    r'(Forward|Reverse)\s+primer\s+\d+\s+\S+\s+\d+\s+Template\s+\d+\s+([.\-A-Za-z]+)\s+\d+',
                    pre):
                key = 'fwd' if a.group(1) == 'Forward' else 'rev'
                aln[key].append(count_mismatch(a.group(2)))
            hits.append({
                'section': section, 'acc': acc, 'desc': desc,
                'gene': extract_gene(desc, acc),
                'product': int(plen.group(1)) if plen else None,
                'fwd_mm': max(aln['fwd']) if aln['fwd'] else None,
                'rev_mm': max(aln['rev']) if aln['rev'] else None,
            })
    return summary_line, hits


def parse_pair_props(html):
    """解析引物对自身的属性表（结果页每对引物一个表）：
    Length、Tm、GC%、Self complementarity、Self 3' complementarity 的 F/R 值。"""
    p = {}
    for key, fld in (('fwd_len', 'FW_PRIMER_LEN'), ('fwd_tm', 'FW_PRIMER_TM'),
                     ('fwd_gc', 'FW_PRIMER_GC'),
                     ('fwd_sc', 'FW_PRIMER_SFCM'), ('fwd_sc3', 'FW_PRIMER_SF3CM'),
                     ('rev_len', 'RV_PRIMER_LEN'), ('rev_tm', 'RV_PRIMER_TM'),
                     ('rev_gc', 'RV_PRIMER_GC'),
                     ('rev_sc', 'RV_PRIMER_SFCM'), ('rev_sc3', 'RV_PRIMER_SF3CM')):
        m = re.search(r'value="([^"]*)"\s+name="%s_\d+"' % fld, html)
        p[key] = m.group(1) if m else None
    return p


def _pair_cells(res):
    """(Length, Tm, GC%, Self complementarity, Self 3' complementarity)
    五格的 'F / R' 展示串。"""
    p = res.get('pair_props') or {}

    def cell(*keys):
        vals = [p.get(k) for k in keys]
        if all(v is None for v in vals):
            return '?'
        return ' / '.join(v if v is not None else '?' for v in vals)
    return (cell('fwd_len', 'rev_len'), cell('fwd_tm', 'rev_tm'),
            cell('fwd_gc', 'rev_gc'), cell('fwd_sc', 'rev_sc'),
            cell('fwd_sc3', 'rev_sc3'))


def _mm(x):
    return 99 if x is None else x


def analyze(hits, target_gene):
    """按基因聚合命中并给出结论。返回 (结论, [(基因, {hits, best}), ...])。"""
    genes = {}
    for h in hits:
        g = genes.setdefault(h['gene'], {'hits': [], 'best': None})
        g['hits'].append(h)
        if g['best'] is None or _mm(h['fwd_mm']) + _mm(h['rev_mm']) < \
                _mm(g['best']['fwd_mm']) + _mm(g['best']['rev_mm']):
            g['best'] = h
    order = sorted(genes.items(),
                   key=lambda kv: (_mm(kv[1]['best']['fwd_mm']) + _mm(kv[1]['best']['rev_mm']),
                                   -len(kv[1]['hits'])))
    target = (target_gene or '').strip()
    if not hits:
        verdict = '无命中'
    elif target:
        off = [g for g in genes if g.lower() != target.lower()]
        if not off:
            verdict = '特异（所有 %d 条产物均在 %s）' % (len(hits), target)
        else:
            verdict = '有脱靶风险：%d 条产物落在 %d 个非目标基因' % (len(hits), len(off))
    elif len(genes) == 1:
        verdict = '单一基因（%s，%d 条产物）' % (next(iter(genes)), len(hits))
    else:
        verdict = '多基因命中：%d 个基因共 %d 条产物（可能存在非特异扩增）' % (len(genes), len(hits))
    return verdict, order


def _name_tags(job):
    """可读文件名片段: (基因或正向引物前6位, 物种)。"""
    gene_tag = re.sub(r'\W+', '_', str(job.get('target_gene') or ''))[:30].strip('_')
    if not gene_tag:
        gene_tag = 'SEQ-' + re.sub(
            r'[^ACGT]', '', (job.get('fwd', '') or '').upper())[:6]
    org_tag = re.sub(r'\W+', '_', str(job.get('organism', 'x')))[:30]
    return gene_tag, org_tag


def run_one(job):
    """提交一对引物做 Primer-BLAST，返回结果字典（失败时含 error 字段）。"""
    res = dict(job)
    res['target_gene'] = (job.get('target_gene') or '').strip()
    try:
        organism = resolve_organism((job.get('organism') or '').strip())
        res['organism_resolved'] = organism
        print('[info] 物种: %s' % organism, file=sys.stderr)
        _, page = submit(build_form(job, organism))
        job_url = find_refresh(page)
        if not job_url:
            err = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', page))
            res['error'] = 'NCBI 未返回任务地址: %s' % err[:300]
            return res
        print('[info] 任务已提交，轮询中...', file=sys.stderr)
        result_url, page = wait_result(job_url, job.get('timeout_min', 10))
        summary, hits = parse_report(page)
        verdict, order = analyze(hits, res['target_gene'])
        res.update(result_url=result_url, summary_line=summary,
                   hits=hits, verdict=verdict, genes=order,
                   pair_props=parse_pair_props(page))
        # head 注入：base 供内嵌/备份时相对路径指回 NCBI；style 隐藏 sandbox 禁 JS
        # 后会常显的帮助说明区（.helpbox/#dnHelp）与 NCBI 页脚
        if re.search(r'<head[^>]*>', page, re.I) and 'skill-hide' not in page:
            inj = ''
            if '<base href' not in page:
                inj += '<base href="https://www.ncbi.nlm.nih.gov/tools/primer-blast/" />'
            inj += ('<style>/* skill-hide */'
                    '.helpbox,#dnHelp,footer{display:none!important}</style>')
            page = re.sub(r'(<head[^>]*>)',
                          lambda m: m.group(1) + '\n' + inj,
                          page, count=1, flags=re.I)
        res['page'] = page  # 内嵌进合并存档报告；默认不再产生单独结果页文件
        if job.get('out'):
            with open(job['out'], 'w', encoding='utf-8') as f:
                f.write(page)
            res['backup'] = os.path.abspath(job['out'])
    except Exception as e:
        res['error'] = '%s: %s' % (type(e).__name__, e)
    return res


def _verdict_class(res):
    if res.get('error'):
        return 'bad'
    v = res.get('verdict', '')
    if v.startswith('无命中'):
        return 'warn'
    if v.startswith('特异') or v.startswith('单一基因'):
        return 'ok'
    return 'bad'


def print_md(res):
    """把单个结果打成给用户看的 Markdown。"""
    print('## Primer-BLAST 引物特异性结果\n')
    print('- 正向引物: `%s` | 反向引物: `%s`' % (
        res.get('fwd', '').strip().upper().replace(' ', ''),
        res.get('rev', '').strip().upper().replace(' ', '')))
    print('- 物种: %s | 数据库: %s' % (
        res.get('organism_resolved', res.get('organism', '?')), res.get('db', 'refseq_mrna')))
    if res.get('error'):
        print('- **错误**: %s\n' % res['error'])
        return
    print('- NCBI 报告状态: %s' % (res.get('summary_line') or '未识别'))
    print('- **原始结果页: %s**' % res.get('result_url', ''))
    if res.get('backup'):
        print('- HTML 备份: %s' % res.get('backup'))
    print('')
    hits = res.get('hits') or []
    if not hits:
        print('### 结论: 无命中\n')
        print('该物种/数据库下，引物对在允许错配范围内没有可扩增产物。'
              '请检查引物方向、物种是否正确，或换用更宽松的数据库（如 core_nt）。')
        return
    print('### 结论: %s\n' % res.get('verdict', ''))
    ln, tm, gc, sc, sc3 = _pair_cells(res)
    print('| 基因 | 产物长度 | F 错配 | R 错配 | Length (F/R) | Tm (F/R) | GC% (F/R) | '
          'Self complementarity (F/R) | '
          "Self 3' complementarity (F/R) | 该基因产物数 |")
    print('|---|---|---|---|---|---|---|---|---|---|')
    target = res.get('target_gene', '')
    for gene, g in res.get('genes') or []:
        b = g['best']
        mark = ''
        if target and gene.lower() != target.lower():
            mark = ' ⚠脱靶'
        elif not target and len(res['genes']) > 1:
            mark = ' ?'
        print('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %d |' % (
            gene + mark, b['product'],
            b['fwd_mm'] if b['fwd_mm'] is not None else '?',
            b['rev_mm'] if b['rev_mm'] is not None else '?',
            ln, tm, gc, sc, sc3, len(g['hits'])))
    print('\n说明: F/R 错配 = 该产物比对中引物的错配碱基数（0 = 完全匹配）；'
          '产物长度取该基因最佳命中（错配最少）产物，多产物条数见末列；'
          "Length/Tm/GC%/Self complementarity/Self 3' complementarity = 引物自身属性（F / R，Self 两项越低越好）；"
          '完整逐条比对见存档报告内嵌结果页。')


def render_html(results, path):
    """生成可长期存档的合并 HTML 报告（原始结果页 URL 会过期，此文件不会）。"""
    css = ('body{font-family:Segoe UI,Arial,sans-serif;max-width:1240px;margin:24px auto;padding:0 16px;color:#222}'
           '.job{border:1px solid #ddd;border-radius:8px;padding:16px;margin:16px 0}'
           '.badge{display:inline-block;padding:2px 10px;border-radius:12px;font-size:13px;font-weight:bold;color:#fff}'
           '.ok{background:#2e9e44}.bad{background:#c93030}.warn{background:#e08a00}'
           'table{border-collapse:collapse;width:100%;margin-top:8px;font-size:14px}'
           'td,th{border:1px solid #ccc;padding:6px 8px;text-align:left}'
           '.meta{color:#666;font-size:13px}.mono{font-family:Consolas,monospace}'
           '.info{color:#333;font-size:16px;line-height:1.7}'
           'details{margin-top:10px}'
           'details summary{cursor:pointer;color:#2b6cb0;font-size:14px;text-decoration:underline}'
           'iframe{width:100%;height:80vh;min-height:640px;border:1px solid #ccc;border-radius:6px;'
           'margin-top:8px;background:#fff}')
    parts = ['<!doctype html><html><head><meta charset="utf-8">'
             '<title>Primer-BLAST 引物特异性存档报告</title><style>%s</style></head><body>' % css,
             '<h1>Primer-BLAST 引物特异性检查报告</h1>',
             '<p class="meta">生成时间: %s | 共 %d 对引物 | 工具: primer-blast-specificity skill'
             '（底层 NCBI Primer-BLAST）</p>' % (time.strftime('%Y-%m-%d %H:%M:%S'), len(results))]
    for i, r in enumerate(results, 1):
        cls = _verdict_class(r)
        title = '错误' if r.get('error') else r.get('verdict', '')
        parts.append('<div class="job"><h2>#%d 引物对 <span class="badge %s">%s</span></h2>' % (
            i, cls, hesc(title)))
        parts.append('<p class="info">正向 <span class="mono">%s</span> / 反向 <span class="mono">%s</span><br>'
                     '物种 %s | 靶基因 %s | 数据库 %s</p>' % (
                         hesc(r.get('fwd', '')), hesc(r.get('rev', '')),
                         hesc(r.get('organism_resolved', r.get('organism', '?'))),
                         hesc(r.get('target_gene') or '未指定'), hesc(r.get('db', 'refseq_mrna'))))
        if r.get('error'):
            parts.append('<p class="meta">错误: %s</p>' % hesc(r['error']))
        if r.get('genes'):
            ln, tm, gc, sc, sc3 = _pair_cells(r)
            parts.append('<table><tr><th>基因</th><th>产物长度</th>'
                         '<th>F 错配</th><th>R 错配</th>'
                         '<th>Length (F/R)</th><th>Tm (F/R)</th><th>GC% (F/R)</th>'
                         "<th>Self complementarity (F/R)</th><th>Self 3' complementarity (F/R)</th>"
                         '<th>该基因产物数</th></tr>')
            target = r.get('target_gene', '')
            for gene, g in r['genes']:
                b = g['best']
                mark = ''
                if target and gene.lower() != target.lower():
                    mark = ' ⚠脱靶'
                elif not target and len(r['genes']) > 1:
                    mark = ' ?'
                parts.append('<tr><td>%s%s</td><td>%s</td><td>%s</td><td>%s</td>'
                             '<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>'
                             '<td>%d</td></tr>' % (
                                 hesc(gene), mark, b['product'],
                                 b['fwd_mm'] if b['fwd_mm'] is not None else '?',
                                 b['rev_mm'] if b['rev_mm'] is not None else '?',
                                 hesc(ln), hesc(tm), hesc(gc), hesc(sc), hesc(sc3),
                                 len(g['hits'])))
            parts.append('</table>')
        if r.get('page'):
            # 完整 NCBI 结果页整页内嵌：单文件存档，URL 过期后仍可离线查看。
            # 存档入口 = 这一个可点开的链接（summary）；临时 URL 不再进报告
            parts.append('<details><summary>点此打开存档的完整 NCBI 结果页</summary>'
                         '<iframe srcdoc="%s" sandbox></iframe></details>' % hesc(r['page']))
        parts.append('</div>')
    parts.append('<p class="meta">F/R 错配 = 比对中引物的错配碱基数（0 = 完全匹配）；'
                 '产物长度取该基因最佳命中（错配最少）产物，多产物条数见末列；'
                 "Length/Tm/GC%/Self complementarity/Self 3' complementarity = 引物自身属性"
                 '（F / R，Self 两项越低越好）。</p>'
                 '</body></html>')
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))


def load_csv(path):
    """读批量 CSV：表头需含 fwd, rev, organism，可选 target_gene（列名不区分大小写）。"""
    jobs = []
    with open(path, newline='', encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            norm = {(k or '').strip().lower(): (v or '').strip() for k, v in row.items()}
            if not norm.get('fwd') or not norm.get('rev') or not norm.get('organism'):
                raise ValueError('CSV 每行必须含 fwd, rev, organism 三列: %s' % row)
            jobs.append({'fwd': norm['fwd'], 'rev': norm['rev'],
                         'organism': norm['organism'],
                         'target_gene': norm.get('target_gene', '')})
    if not jobs:
        raise ValueError('CSV 没有数据行')
    return jobs


def main():
    p = argparse.ArgumentParser(description='NCBI Primer-BLAST 引物特异性检查')
    p.add_argument('--fwd', default='', help="正向引物 5'->3'（单对模式必填）")
    p.add_argument('--rev', default='', help="反向引物 5'->3'（单对模式必填）")
    p.add_argument('--organism', default='', help='物种名或含 taxid 的字符串（缺省按人源 Homo sapiens）')
    p.add_argument('--csv', default='', help='批量模式: CSV 文件，表头 fwd,rev,organism[,target_gene]')
    p.add_argument('--report', default='', help='存档 HTML 报告路径（单对/批量通用，默认按时间命名）')
    p.add_argument('--delay', type=float, default=5, help='批量模式每对引物提交间隔秒数')
    p.add_argument('--db', default='refseq_mrna',
                   help='refseq_mrna(默认)/refseq_representative_genomes/core_nt/refseq_rna/nt')
    p.add_argument('--target-gene', default='', help='预期靶基因符号（可选），用于判定脱靶')
    p.add_argument('--product-min', type=int, default=70)
    p.add_argument('--product-max', type=int, default=1000)
    p.add_argument('--mismatch', type=int, default=2, help='允许的最大总错配（NCBI 最低 2）')
    p.add_argument('--mismatch-3end', type=int, default=2, help='3\' 端允许最大错配')
    p.add_argument('--num-return', type=int, default=10)
    p.add_argument('--entrez-query', default='')
    p.add_argument('--out', default='',
                   help='可选: 额外保存一份单独的结果页 HTML（默认不保存，完整结果页已内嵌存档报告）')
    p.add_argument('--timeout-min', type=int, default=10)
    args = p.parse_args()

    common = dict(db=args.db, product_min=args.product_min, product_max=args.product_max,
                  mismatch=args.mismatch, mismatch_3end=args.mismatch_3end,
                  num_return=args.num_return, entrez_query=args.entrez_query,
                  timeout_min=args.timeout_min)

    if args.csv:
        try:
            jobs = load_csv(args.csv)
        except Exception as e:
            print('CSV 读取失败: %s' % e, file=sys.stderr)
            sys.exit(2)
        print('[info] 批量任务 %d 对，提交间隔 %gs' % (len(jobs), args.delay), file=sys.stderr)
        results = []
        for i, j in enumerate(jobs):
            print('\n--- [%d/%d] ---' % (i + 1, len(jobs)))
            res = run_one(dict(common, **j))
            results.append(res)
            print_md(res)
            if i < len(jobs) - 1:
                time.sleep(args.delay)
        report = args.report or os.path.abspath(
            'primer_blast_report_%s.html' % time.strftime('%Y%m%d_%H%M%S'))
        render_html(results, report)
        nerr = sum(1 for r in results if r.get('error'))
        print('\n批量完成: %d 成功 / %d 失败（共 %d）' % (len(results) - nerr, nerr, len(results)))
        print('存档报告: %s' % os.path.abspath(report))
        sys.exit(1 if nerr == len(results) else 0)

    if not (args.fwd and args.rev):
        p.error('单对模式需 --fwd/--rev；批量模式用 --csv')
    organism = args.organism or 'Homo sapiens'   # 未给种属默认人源
    if not args.organism:
        print('[info] 未指定 --organism，按默认人源 Homo sapiens 运行', file=sys.stderr)
    res = run_one(dict(common, fwd=args.fwd, rev=args.rev, organism=organism,
                       target_gene=args.target_gene, out=args.out))
    print_md(res)
    if not res.get('error'):
        # 单对模式同样生成带徽章的存档报告
        gt, ot = _name_tags({'target_gene': args.target_gene,
                             'fwd': args.fwd, 'organism': organism})
        report = args.report or os.path.abspath(
            'PrimerBlast_report_%s_%s_%s.html' % (
                gt, ot, time.strftime('%Y%m%d_%H%M%S')))
        render_html([res], report)
        print('存档报告: %s' % os.path.abspath(report))
    sys.exit(2 if res.get('error') else 0)


if __name__ == '__main__':
    main()
