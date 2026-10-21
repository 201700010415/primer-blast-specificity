# Primer-BLAST Specificity Check

**给一对引物，还你一个敢写进论文的结论。**

提交 NCBI Primer-BLAST → 自动轮询等结果 → 用红/黄/绿结论告诉你这对引物是**特异**、**多基因命中**还是**有脱靶风险**，并生成可离线存档的 HTML 报告。顺手还能按基因名从 PrimerBank 查现成的 qPCR 引物，一键串联验证。

![Python 3](https://img.shields.io/badge/Python-3-blue)
![No dependencies](https://img.shields.io/badge/dependencies-none-brightgreen)
![NCBI Primer-BLAST](https://img.shields.io/badge/NCBI-Primer--BLAST-326CE5)
![PrimerBank](https://img.shields.io/badge/PrimerBank-qPCR-e6522c)

> 本仓库同时是一个 Agent Skill：`SKILL.md` 定义了完整的使用与结果解读规范，可直接接入 AI 编程助手自动执行。

## 为什么需要它

拿到一对引物，最怕的就是非特异扩增——跑完实验才发现 Cq 曲线不对，几周就没了。

手动在 NCBI 网页上操作 Primer-BLAST：填表、提交、等 1–2 分钟、逐条翻命中基因……一对引物还好，十几对引物就是半天。这个工具把整条链路做成了命令：

- **单对**：一条命令出结论
- **批量**：一个 CSV 全部提交，合并成一份带红黄绿徽章的存档报告
- **查引物**：按基因名（支持别名）直接从 PrimerBank 拿经过验证的 qPCR 引物对

## 特性

- **一眼看懂的结论** —— 绿=特异/单一基因，黄=多基因命中，红=有脱靶风险，不用自己翻比对表
- **逐基因错配明细** —— 每个命中基因一行：产物长度、正/反向引物错配数、Tm、GC%、自互补性
- **报告可离线存档** —— 生成的 HTML 内嵌完整 NCBI 结果页，临时链接过期后仍能逐位点比对
- **PrimerBank 引物查询** —— 本地内置人+鼠 30.6 万条基因名/GeneID 映射（NCBI gene_info 权威数据），官方 symbol、别名、GeneID 都能查
- **一键串联** —— `--blast` 查到引物后自动逐对跑特异性验证
- **零依赖** —— 纯 Python 3 标准库，克隆下来就能跑，不用 pip 装任何东西

## 快速开始

需要 Python 3 和对 `ncbi.nlm.nih.gov` 的网络访问，仅此而已。

### 1. 单对引物

```bash
python scripts/blast_primer.py \
  --fwd ATGCCACCTGCTGAACGGATT \
  --rev CTCCGGAAGAAACTGCACCA \
  --organism "Homo sapiens" \
  --target-gene HMOX1
```

`--organism` 省略时默认人源 `Homo sapiens`；给了 `--target-gene` 才能明确判定"脱靶"，否则按命中基因数提示。

### 2. 批量验证

准备一个 CSV（表头必需，列名不区分大小写）：

```csv
fwd,rev,organism,target_gene
ATGCCACCTGCTGAACGGATT,CTCCGGAAGAAACTGCACCA,Homo sapiens,HMOX1
ACAGGCTATGGCTTCCAGGT,MCTTCTGCACCTTCGACCAT,Mus musculus,Hmox1
```

```bash
python scripts/blast_primer.py --csv primers.csv --report my_report.html
```

多对引物、多物种一律用批量模式，脚本会自动控制提交间隔、失败重试，最后汇总。

### 3. 查 PrimerBank 现成引物

```bash
python scripts/primerbank_query.py --gene HMOX1 --species mouse
```

结果自动保存为 `HMOX1_Mus_musculus.md`（stdout 同步输出）。`--gene` 支持官方 symbol、**别名**、纯数字 GeneID——比如查 `HMOX1` 的鼠源引物，输入别名 `HO1` 也能解析。

### 4. 查完直接验（一条龙）

```bash
python scripts/primerbank_query.py --gene HMOX1 --species mouse --blast
```

拿到引物对后自动逐对提交 Primer-BLAST，出每对结论 + 合并存档报告 `PrimerBlast_report_HMOX1_Mus_musculus_<时间>.html`。

## 结果长什么样

**终端输出（stdout，Markdown）**：结论行 + 每基因一张表

| 产物长度 | F/R 错配 | Tm (F/R) | GC% (F/R) | Self complementarity (F/R) | 产物数 |
|---|---|---|---|---|---|
| 152 bp | 0 / 1 | 60.1 / 59.8 | 55 / 52 | 4.0 / 3.0 | 1 |

- **F/R 错配 = 该产物比对中引物的错配碱基数**，0 为完全匹配，错配越多结合越弱
- 产物长度只显示该基因最佳命中（错配最少）那条；同基因多产物看末列
- Self complementarity 两项越低越好（引物自身二聚体风险）

**存档 HTML 报告**（唯一需要留档的结果文件）：

```
PrimerBlast_report_<基因|引物前6位>_<物种>_<时间>.html
```

每对引物一张卡片：结论徽章 → 引物/物种信息 → 统计表 → 可折叠的"完整 NCBI 结果页存档"（点开就地展开，离线可看，联网可恢复官方样式）。

> 原始结果页 URL 只出现在终端输出里，不进存档报告——临时链接会过期，**存档报告才是最终交付物**。

**PrimerBank 查询记录**：`<基因名>_<物种>.md`，重复查询自动覆盖刷新。

## 主要参数（blast_primer.py）

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--fwd` / `--rev` | — | 正/反向引物（5'→3'，仅碱基），单对模式必填 |
| `--organism` | `Homo sapiens` | 物种名，支持 `taxid:9606` 写法 |
| `--csv` | — | 批量模式 CSV 路径，表头 `fwd,rev,organism[,target_gene]` |
| `--target-gene` | — | 预期靶基因，用于判定脱靶 |
| `--db` | `refseq_mrna` | 基因组水平换 `refseq_representative_genomes`，跨物种/未知换 `core_nt` |
| `--report` | 按时间命名 | 存档 HTML 报告路径 |
| `--mismatch` / `--mismatch-3end` | `2` / `2` | 允许的最大总错配 / 3' 端错配（NCBI 最低 2） |
| `--product-min` / `--product-max` | `70` / `1000` | 产物长度范围（bp） |
| `--delay` | `5` | 批量模式每对提交间隔（秒） |
| `--timeout-min` | `10` | 单任务轮询超时（分钟） |

进度日志走 stderr，可忽略；stdout 是给人看的 Markdown 报告，可直接重定向存档。

## 数据来源与合规

- `data/gene_map.tsv` 来自 NCBI **gene_info**（人+鼠，含 Synonyms 列），查询超过 180 天脚本会自动提醒重跑 `scripts/build_gene_map.py`（约 9MB 下载）
- 引物结果不缓存，每次实时查询 PrimerBank
- NCBI 请求已携带 `tool` 标识；设置环境变量 `NCBI_EMAIL` 可自动附带 `email` 字段，更符合 NCBI ToU 要求
- 单个任务约 30–120 秒，脚本内部自动轮询与重试（瞬时网络错误重试 3 次）

## 局限

- **只做特异性检查，不做引物设计**（那是 Primer-BLAST 的设计模式，本脚本固定用你给定的引物）
- NCBI 页面改版可能导致解析失败；届时把结果页 HTML 存下来对照修复 `scripts/blast_primer.py` 的 `parse_report`
- PrimerBank 覆盖以人/鼠为主，未收录的引物会如实说明
- 基因名匹配到多个 GeneID 或无匹配时会列出候选让你选，不会替你猜

## 目录结构

```
primer-blast-specificity/
├── SKILL.md                      # Agent Skill 定义（使用规范 + 结果解读铁律）
├── README.md
├── scripts/
│   ├── blast_primer.py           # Primer-BLAST 提交 / 轮询 / 解析 / 报告生成
│   ├── primerbank_query.py       # PrimerBank 基因引物查询（--blast 串联验证）
│   └── build_gene_map.py         # 重建本地基因名→GeneID 映射表
└── data/
    └── gene_map.tsv              # 人+鼠 30.6 万条基因映射（NCBI gene_info）
```
