# GTDB Renew

`GTDB Renew` 现在分成三部分，统一通过 `gtdbkit` 管理：

1. **网页版**
   - GTDB 版本比较
   - 网页可视化浏览
   - 网页下载与比较
2. **命令行下载工具**
   - 按 GTDB `taxon` 下载
   - 按 `accession` 下载
   - 使用 GTDB metadata 提取 accession
   - 使用 `ncbi_datasets` 分批并行下载
3. **数据库构建与维护工具**
   - 基因组文件重命名
   - contig / sequence header 统一
   - CheckM2 质控
   - barrnap 提取 rRNA
   - 16S 等 rRNA 序列输出
   - 版本化数据库管理

## 新入口名称

统一入口现在推荐使用：

```bash
pixi run gtdbkit
```

也提供一个更短的别名：

```bash
pixi run gtk
```

旧入口 `gttdb` 仍然保留兼容，但以后建议都使用 `gtdbkit`。

## 平台支持

| 平台 | 状态 |
|---|---|
| macOS arm64 | ✅ 已测试 |
| Linux x86_64 | ✅ 已测试 |
| Linux aarch64 | ⏳ 理论可用 |
| Windows | ❌ 不支持 |

## 环境安装

项目使用 [Pixi](https://pixi.sh) 管理环境和依赖。

```bash
git clone https://github.com/Changhai996/GTDB_genome_download
cd GTDB_genome_download
pixi install
```

如果是 Linux / macOS，也可以直接：

```bash
chmod +x run.sh
./run.sh help
```

## 第一部分：网页版

打开网页版：

```bash
pixi run gtdbkit web
```

或：

```bash
pixi run web
```

网页版主文件：

- `app.py`

主要功能：

- GTDB 版本比较
- Custom Download
- Dataset Updater
- Database Builder

## 第二部分：GTDB 下载

主命令：

```bash
pixi run gtdbkit fetch --help
```

### 1. 按 taxon 下载

```bash
pixi run gtdbkit fetch \
  -R 220.0 \
  -t "c__Bathyarchaeia" \
  -m representative \
  -j 4 \
  -b 50
```

说明：

- `-t / --taxon`：指定 GTDB 类群
- `-m / --scope`：`representative`、`all`、`accessions-only`
- `-j / --threads`：并行下载批次数
- `-b / --batch-size`：每个 `datasets` 批次的 accession 数量

### 2. 按 accession 下载

```bash
pixi run gtdbkit fetch \
  -a GCA_003151735.1,GCF_000196175.1 \
  -j 4 \
  -b 20
```

或：

```bash
pixi run gtdbkit fetch -A accession_list.txt
```

### 3. 直接导入第三部分输入目录

第二部分下载的数据，可以直接作为第三部分输入：

```bash
pixi run gtdbkit fetch \
  -t "c__Bathyarchaeia" \
  -i /Users/duanchanghai/Downloads/Database/GTDB_Bathyarchaeia
```

`-i / --import-dir` 支持两种方式：

- 默认 `symlink`
- 可选 `copy`

## 第三部分：数据库构建与维护

主命令：

```bash
pixi run gtdbkit build-db --help
```

### 1. 以总目录方式构建

```bash
pixi run gtdbkit build-db \
  -n Bathyarchaeia \
  -v v1.0 \
  -D /Users/duanchanghai/Downloads/Database \
  -o local_databases \
  -j 8 \
  -Q \
  -B \
  -c /path/to/uniref100.KO.1.dmnd
```

说明：

- `-D / --source-root`：总输入目录，一级子目录视为来源
- `-n / --name`：数据库名称
- `-v / --version`：数据库版本
- `-Q / --checkm2`：启用 CheckM2
- `-B / --barrnap`：启用 barrnap
- `-c / --checkm2-db`：CheckM2 数据库路径
- `--exclude-hidden / --no-exclude-hidden`：是否跳过 `.DS_Store`、`._*` 等隐藏文件
- `--strict-fasta-check / --no-strict-fasta-check`：开启时按文件内容校验 FASTA，并自动识别无标准后缀的 FASTA；关闭时回退为仅按后缀扫描
- `-S / --collect-16s-to`：把数据库里所有 16S 汇总到一个输出 FASTA

### 2. 以多个来源目录方式构建

```bash
pixi run gtdbkit build-db \
  -n Bathyarchaeia \
  -v v1.0 \
  -s /path/to/source_A /path/to/source_B \
  -j 8 \
  -Q \
  -B
```

### 第三部分会做什么

- 递归扫描 `fa / fna / fasta`，默认还会自动识别无标准后缀但内容符合 FASTA 格式的文件
- 不修改原始文件
- 标准化文件名
- 标准化内部序列 header
- contig header 与最终 genome 名统一
- 输出重命名前后映射表
- 统计 `Genome_Size_bp / Contig_Count / GC_Content_%`
- 运行 CheckM2，输出 `Completeness / Contamination / GTDB_Score`
- 运行 barrnap，输出 GFF 与 rRNA fasta
- 额外提取 16S rRNA fasta，header 保留修改后的基因组名，并在需要时追加序号
- 可把数据库中全部 16S 再汇总为一个总 FASTA 文件
- 生成 metadata、log、版本比较表
- 支持基于 MD5 的增量版本处理
- 支持按步骤断点续传：若当前版本的标准化基因组、barrnap、16S 或 CheckM2 结果已存在，会优先复用
- 启动时打印“隐藏/异常 FASTA 共跳过多少个”以及隐藏文件、异常 FASTA、无后缀 FASTA 的发现统计

### 3. 仅汇总当前数据库全部 16S

如果当前版本数据库已经构建完成，只想把全部 16S 提取到一个总文件，可以直接：

```bash
pixi run gtdbkit build-db \
  -n Bathyarchaeia \
  -v v1.0 \
  -o local_databases \
  -S /path/to/Bathyarchaeia_all_16S.fasta
```

说明：

- 这种用法不需要再次提供 `-D` 或 `-s`
- 程序会直接复用当前版本目录下已有的 `metadata`、标准化基因组和 `barrnap_results`
- 如果单个基因组缺少 `*_16S.fasta`，但已有对应 GFF，则会先从现有 GFF 补提取，再汇总

## 第二 + 第三部分一体化

如果想“先下载 GTDB，再自动导入并构库”，使用：

```bash
pixi run gtdbkit prepare-db \
  -n Bathyarchaeia \
  -v v1.0 \
  -t "c__Bathyarchaeia" \
  -i /Users/duanchanghai/Downloads/Database/GTDB_Bathyarchaeia \
  -j 8 \
  -Q \
  -B \
  -c /path/to/uniref100.KO.1.dmnd
```

注意：

- 这里需要 `-t` 或 `-a`，是因为第一步要先下载
- 并不是第三部分本身需要 `taxon`

## CheckM2 数据库路径

构建阶段会按以下优先级寻找 CheckM2 DIAMOND 数据库：

1. `-c / --checkm2-db`
2. 环境变量 `CHECKM2DB` / `DIAMOND_DB`
3. `<project>/checkm2_database/CheckM2_database/*.dmnd`
4. `<project>/checkm2_database/*.dmnd`

下载 CheckM2 数据库：

```bash
pixi run checkm2 database --download
```

## 典型输出目录

```text
local_databases/Bathyarchaeia_v1.0/
├── genomes/
├── barrnap_results/
├── checkm2_results/
├── Bathyarchaeia_v1.0_metadata.csv
├── Bathyarchaeia_v1.0_genome_id_mapping.csv
├── Bathyarchaeia_v1.0_scan_inventory.csv
├── Bathyarchaeia_v1.0_source_counts.csv
├── Bathyarchaeia_v1.0_version_comparison.csv
└── build_summary.log
```

## 主要文件

```text
GTDB_renew/
├── app.py
├── gtdbkit.py
├── gttdb.py
├── gtdb_downloader.py
├── db_builder_cli.py
├── checkm2_compat/
│   └── sitecustomize.py
├── pixi.toml
├── run.sh
└── README.md
```

说明：

- `gtdbkit.py`：新的统一入口
- `gttdb.py`：旧入口兼容包装器
- `gtdb_downloader.py`：第二部分下载器
- `db_builder_cli.py`：第三部分数据库构建器

## 兼容旧命令

这些旧名字仍可用，但不再推荐：

- `gttdb` -> `gtdbkit`
- `download-gtdb` -> `fetch`
- `build` -> `build-db`
- `run-all` -> `prepare-db`

## License

内部工具；`CheckM2`、`barrnap`、`DIAMOND` 遵循各自许可证。
