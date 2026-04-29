# GTDB Genome Download & Taxonomy Analyzer

easily download genome from each GTDB release

This tool is a comprehensive Python-based toolkit designed for researchers working with the **Genome Taxonomy Database (GTDB)** and **NCBI Genomes**. It provides an intuitive Web GUI (built with Streamlit) to automatically download GTDB releases, analyze taxonomic changes across versions, and bulk-download genome sequences (`.fna`) directly from NCBI with a robust resume capability.

## Key Features

1. **Automated GTDB Data Retrieval**
   - Automatically crawls and downloads all historical and latest Archaea taxonomy files (`ar*_taxonomy_r*.tsv`) from the official GTDB releases server.

2. **Single Version Explorer**
   - Browse the taxonomy structure of any specific GTDB release.
   - Search for specific taxonomic groups (e.g., `c__Bathyarchaeia`) and instantly view all associated genomes.

3. **Version Comparison & Tracking**
   - Compare two GTDB releases to identify newly added, removed, or reclassified genomes within a specific clade.
   - Generates detailed statistics on taxonomic units (Phylum, Class, Order, Family, Genus, Species) and net genome count changes.

4. **Robust NCBI Datasets Integration**
   - Automatically installs the correct `ncbi_datasets` CLI tool for your OS (macOS/Linux).
   - **Batch Downloading & Resume:** Downloads genomes in small batches (e.g., 20 at a time) to prevent large zip failures. Automatically skips already downloaded files (断点续传) so you never lose progress if the network drops.
   - Supports one-click downloading of "All genomes in a group" or "Only newly added genomes" from the comparison tab.

5. **Custom Download & Email Reporting**
   - Paste a custom list of Genome IDs (with or without `RS_`/`GB_` prefixes) to download their `.fna` sequences.
   - Automatically generates a `.csv` report detailing successful and failed downloads.
   - Configurable SMTP settings to automatically send the generated report as an email attachment to your inbox (e.g., `changhaiduan@163.com`).

---

## Installation

### Prerequisites
- Python 3.8 or higher
- Git

### Setup
1. Clone this repository:
   ```bash
   git clone https://github.com/Changhai996/GTDB_genome_download.git
   cd GTDB_genome_download
   ```

2. Create a virtual environment and activate it:
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate
   ```

3. Install the required dependencies:
   ```bash
   pip install pandas streamlit requests beautifulsoup4
   ```

---

## Usage

### 1. Download GTDB Taxonomy Data
Before running the analyzer, you need to download the raw taxonomy files from GTDB. Run the included script:

```bash
python download_gtdb_taxonomy.py
```
This script will create a `gtdb_data` folder and download `.tsv` files from the GTDB releases server.

### 2. Launch the Web Application
Start the Streamlit application:

```bash
streamlit run app.py
```
This will open the application in your default web browser (usually at `http://localhost:8501`).

### 3. Navigating the App
The application is divided into three main tabs:
- **🔍 Single Version Explorer:** Select a version, search for a taxon, and download all its genomes.
- **📊 Version Comparison:** Select multiple versions, analyze the differences in a specific taxon, and download either the full set or just the newly added genomes.
- **📥 Custom Download & Email:** Paste custom Genome IDs, configure your SMTP email settings, and receive a detailed download status report straight to your inbox.

---

## File Structure

- `app.py`: The main Streamlit web application and NCBI download logic.
- `download_gtdb_taxonomy.py`: The scraper script to fetch GTDB taxonomy TSV files.
- `gtdb_data/`: Directory where GTDB taxonomy files are stored.
- `ncbi_downloads/`: Directory where downloaded `.fna` genome sequences and reports are saved.
- `datasets`: The NCBI datasets CLI executable (automatically downloaded by `app.py`).

## Notes
- The NCBI `datasets` tool downloads genome sequences. If some genomes fail to download, it may be because they have been suppressed or removed from NCBI. The tool will generate a report listing these failed IDs.
- For the email feature, ensure your email provider allows SMTP connections. For services like 163 or Gmail, you usually need to generate an "Authorization Code" or "App Password" instead of using your primary login password.
