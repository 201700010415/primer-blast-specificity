# -*- coding: utf-8 -*-
"""构建基因名 <-> GeneID 本地映射表（人 + 鼠）。

数据源: NCBI FTP 分物种 gene_info（权威、含别名）。
输出:   ../data/gene_map.tsv  (tax_id, gene_id, symbol, synonyms, description)
重跑即刷新；primerbank_query.py 查询用它把基因名解析成 GeneID。
"""
import csv
import gzip
import io
import os
import sys
import urllib.request

SOURCES = [
    ('9606', 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'),
    ('10090', 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz'),
]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'gene_map.tsv')


def fetch_gz(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'primer-blast-specificity skill'})
    data = urllib.request.urlopen(req, timeout=300).read()
    return gzip.GzipFile(fileobj=io.BytesIO(data)).read().decode('utf-8', 'ignore')


def main():
    rows = []
    for tax_id, url in SOURCES:
        print('[info] 下载 %s' % url.rsplit('/', 1)[-1], file=sys.stderr)
        text = fetch_gz(url)
        reader = csv.reader(io.StringIO(text), delimiter='\t')
        header = None
        n = 0
        for rec in reader:
            if not rec:
                continue
            if rec[0].startswith('#'):
                header = [c.lstrip('#') for c in rec]
                col = {name: i for i, name in enumerate(header)}
                continue
            gene_id = rec[col['GeneID']]
            symbol = rec[col['Symbol']]
            syn = rec[col['Synonyms']] if 'Synonyms' in col else '-'
            desc = rec[col['description']] if 'description' in col else ''
            rows.append((tax_id, gene_id, symbol, syn, desc))
            n += 1
        print('[info] taxid %s: %d 条' % (tax_id, n), file=sys.stderr)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['tax_id', 'gene_id', 'symbol', 'synonyms', 'description'])
        w.writerows(rows)
    print('[info] 映射表已写入: %s（共 %d 条）' % (os.path.abspath(OUT), len(rows)), file=sys.stderr)


if __name__ == '__main__':
    main()
