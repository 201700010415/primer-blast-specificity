---
name: primer-blast-specificity
description: Submit a primer pair plus organism to NCBI Primer-BLAST, poll for the report, and summarize per-gene hits with forward/reverse mismatch counts plus the raw result URL. Also fetches ready-made qPCR primer pairs from PrimerBank by gene (local gene name/synonym to GeneID map). Use when the user asks to verify primer specificity, check primer off-target binding, 验证/检验引物特异性, 引物脱靶检查, Primer-BLAST a primer pair, or 查 PrimerBank 某基因的引物 (e.g. "查 HMOX1 鼠源引物").
---

# Primer-BLAST Specificity Check

给定一对引物（5'->3'）和物种名，调用 NCBI Primer-BLAST 检查特异性，输出：结论、每个基因的最佳命中与 F/R 错配数、原始结果页 URL、本地 HTML 备份。

## 输入要求

- 必需：正向引物、反向引物（仅碱基序列）
- 可选：物种名（如 `Homo sapiens`，或带 `taxid:9606`）——**未给时默认人源 `Homo sapiens`**，脚本会提示默认值
- 可选：预期靶基因 `--target-gene`（有它才能明确判定"脱靶"；没有则按"命中基因数"判断）
- 可选：数据库 `--db`（默认 `refseq_mrna`；基因组水平可换 `refseq_representative_genomes`，跨物种/未知可换 `core_nt`）

## 运行

单对：

```bash
python "<skill目录>/scripts/blast_primer.py" --fwd <正向引物> --rev <反向引物> --organism "<物种>" [--target-gene <基因>] [--db <数据库>]
```

批量（多对引物或多物种一律用批量模式，不要逐个手动循环）：

```bash
python "<skill目录>/scripts/blast_primer.py" --csv <CSV文件> [--report <报告.html>] [--delay 5]
```

- CSV 表头必需 `fwd,rev,organism`，可选 `target_gene`（列名不区分大小写，UTF-8）。缺表头或列时帮用户生成 CSV 再运行。
- 批量模式额外产出一份**合并存档 HTML 报告**（结论带红/黄/绿徽章 + 各基因错配表 + 结果页链接），`--report` 指定路径，默认写到当前目录按时间命名。
- 仅用 Python 3 标准库；每个任务约 30-120 秒，脚本内部已轮询，`--timeout-min` 默认 10。
- 脚本 stdout 是给用户看的 Markdown 报告，stderr 是进度日志，可忽略。

### 等待阶段的汇报节奏（重要）

BLAST 提交后的等待阶段**不要直接结束回合**，要 hold 在当前回合内主动盯：

- **每 30 秒查一次**任务状态（进程是否结束、报告文件是否生成），用短命令检查，不要傻等，从提交起主动盯到 **10 min**（与 `--timeout-min` 默认上限一致）。
- 查到结果已出：**不要停下问用户，直接往下走**——读取输出、生成存档报告，汇报结论 + 报告链接。
- 超过 **10 min** 仍无结果才停：结束本轮等待，但**不要杀任务**，让它在后台继续跑；收到完成通知后同样直接继续出报告并汇报。

## PrimerBank 基因引物查询（按基因名/GeneID 拿现成 qPCR 引物）

```bash
python "<skill目录>/scripts/primerbank_query.py" --gene <基因名或GeneID> --species human|mouse [--csv out.csv] [--blast]
```

- 输入支持：官方基因 symbol、**别名**、纯数字 GeneID。解析靠本地映射表 `data/gene_map.tsv`（人+鼠 30.6 万条，NCBI gene_info 权威数据，含 Synonyms）。
- 报告**自动保存本地**：文件名 = `基因名_种属.md`（如 `Hmox1_Mus_musculus.md`），写到运行时的当前目录；stdout 同步输出同样内容。重复查询同名文件直接覆盖刷新。
- 查询链路：基因名 --(本地映射表)--> GeneID --(官方 symbol 提交)--> PrimerBank 在线接口。注：PrimerBank 的 "NCBI Gene ID" 搜索方式服务端已失效（2026-10 实测返回空页），故用官方 symbol 提交；结果页显示的 GeneID 与本地解析自动核对，不一致会告警。
- 引物结果**不缓存**（每次实时查）；映射表建一次反复用，超过 180 天脚本自动提醒重跑 `scripts/build_gene_map.py`（约 9MB 下载）。
- `--csv` 导出 `fwd,rev,organism,target_gene` 格式，正好是下方批量 blast 的输入。
- **`--blast` 一键串联**：查到引物后自动逐对跑 Primer-BLAST 特异性（间隔 `--delay` 秒，默认 5），出每对结论 + 合并存档报告 `PrimerBlast_report_<基因>_<物种>_<时间>.html`。用户只要求查询时不加该参数；用户要"验证/查特异性"时加。预估耗时：每对约 30-120s。
- 基因名匹配到多个 GeneID 或无匹配时：列出候选/提示拼写，让用户选定后再查，不要自己猜。
- 引物未收录时如实说明（PrimerBank 覆盖人/鼠为主）。

## 结果解读与呈现

**如实反馈（铁律）**：只汇报脚本实际解析到的结果——解析失败、无命中、超时就如实说，绝不编造、不凭猜测补数据、不擅自换基因/物种重跑凑结果；遇到歧义（如基因匹配多个 GeneID）停下来列出候选问用户，用户没确认前不要自己猜着往下跑。

1. **结论行**：
   - `特异/单一基因`：所有产物在同一基因。
   - `有脱靶风险`：产物落在目标基因之外的基因，列出的 ⚠脱靶 行要醒目提示。
   - `多基因命中`（未给 target-gene）：提示可能非特异扩增，建议用户提供预期靶基因复跑。
   - `无命中`：引物在该物种/数据库下无可扩增产物，提示检查引物方向、物种拼写或换数据库。
2. **表格**（Markdown 与存档 HTML 一致）：每基因一行，列为 产物长度、F/R 错配、Length (F/R)、Tm (F/R)、GC% (F/R)、Self complementarity (F/R)、Self 3' complementarity (F/R)、该基因产物数（不列 accession，逐位点比对看内嵌结果页）。F/R 错配 = 该产物比对中引物的错配碱基数，0 = 完全匹配，错配越多结合越弱；**产物长度只显示该基因最佳命中（错配最少）那条产物**，同基因多产物的条数看末列；后 5 列是引物自身属性（F / R 并列，其中两个 Self complementarity 越低越好）。原始结果页 URL 只出现在终端 Markdown 输出里，不进存档报告（临时链接会过期，最终交付物是存档报告本身）。
3. **存档报告是唯一结果文件**（单对与批量均自动生成）：`PrimerBlast_report_<基因|引物前6位>_<物种>_<时间>.html`（绿=特异/单一基因，黄=多基因命中，红=有脱靶风险）。每对引物卡片自上而下：结论徽章 → 引物/物种信息段（16px 大字，不显示 NCBI 状态行）→ **统计表** → 存档入口。存档入口为唯一的链接式折叠区"点此打开存档的完整 NCBI 结果页"（在统计表下方，折叠区内只有内嵌页，无其他链接），点开即就地展开内嵌存档页（存档不过期，URL 失效后仍可离线查看逐位点比对，联网打开可恢复官方样式）。内嵌页已自动隐藏 NCBI 帮助说明区和页脚。**不再产生单独的结果页备份文件**；确需单独文件时才加 `--out <路径>`。

## 限制与排错

- NCBI 页面改版可能导致解析失灵：脚本失败时把提交页/结果页 HTML 保存下来，检查 `prPairTl/prPairDtl` 结构是否变化，修复 `scripts/blast_primer.py` 中 `parse_report` 的正则前先与用户确认。
- 提交/轮询失败：脚本已对瞬时网络错误（NCBI 偶发超时、限速）自动重试 3 次；仍失败则该条记为错误并继续批量，最后汇总失败数。
- NCBI 请求已自动携带 `tool` 标识；如需符合 NCBI ToU 的完整标识，设置环境变量 `NCBI_EMAIL`（设了会自动附带 `email` 字段）。
- 本 skill 只做引物对特异性检查；不做引物设计（那是 Primer-BLAST 的设计模式，本脚本固定用给定引物）。
