import os
import re
import pandas as pd
import streamlit as st
import platform
import stat
import urllib.request
import subprocess
import zipfile
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

st.set_page_config(page_title="GTDB Taxonomy Analyzer", layout="wide")

@st.cache_data
def load_data(data_dir="gtdb_data"):
    all_data = []
    if not os.path.exists(data_dir):
        return pd.DataFrame()
        
    files = os.listdir(data_dir)
    tsv_files = [f for f in files if f.endswith('.tsv')]
    
    for f in tsv_files:
        match = re.search(r'_r(\d+)\.tsv', f)
        if not match:
            continue
        version = int(match.group(1))
        filepath = os.path.join(data_dir, f)
        df = pd.read_csv(filepath, sep='\t', header=None, names=['Genome_ID', 'Taxonomy'])
        df['Version'] = version
        tax_levels = df['Taxonomy'].str.split(';', expand=True)
        for i, level in enumerate(['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']):
            if i < tax_levels.shape[1]:
                df[level] = tax_levels[i]
            else:
                df[level] = None
        all_data.append(df)
        
    if not all_data:
        return pd.DataFrame()
        
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df

def ensure_datasets_cli():
    if os.path.exists('datasets'):
        return True
    system = platform.system().lower()
    if system == 'darwin':
        url = 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
    elif system == 'linux':
        url = 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
    else:
        return False
    try:
        urllib.request.urlretrieve(url, 'datasets')
        st_mode = os.stat('datasets').st_mode
        os.chmod('datasets', st_mode | stat.S_IEXEC)
        return True
    except Exception as e:
        st.error(f"Failed to download NCBI datasets CLI: {e}")
        return False

def get_downloaded_accessions(output_dir):
    downloaded = set()
    if not os.path.exists(output_dir):
        return downloaded
    for file in os.listdir(output_dir):
        if file.endswith('.fna'):
            match = re.match(r'^(GCA_\d+\.\d+|GCF_\d+\.\d+)', file)
            if match:
                downloaded.add(match.group(1))
    return downloaded

def download_genomes(genome_list, output_dir, zip_name="dataset.zip", progress_bar=None, status_text=None, batch_size=20):
    if not genome_list:
        return False, "No genomes to download.", [], []
        
    all_accessions = [g.strip()[3:] if g.strip().startswith(("RS_", "GB_")) else g.strip() for g in genome_list if g.strip()]
    total_genomes = len(all_accessions)
    
    os.makedirs(output_dir, exist_ok=True)
    
    downloaded_accs = get_downloaded_accessions(output_dir)
    pending_accessions = [acc for acc in all_accessions if acc not in downloaded_accs]
    
    if status_text:
        if downloaded_accs:
            status_text.info(f"Found {len(downloaded_accs)} already downloaded genomes. {len(pending_accessions)} remaining out of {total_genomes}.")
        else:
            status_text.info(f"Preparing to download {total_genomes} genomes into `{output_dir}` in batches of {batch_size}...")
            
    if not pending_accessions:
        if progress_bar:
            progress_bar.progress(1.0)
        return True, f"All {total_genomes} genomes have already been downloaded to `{os.path.abspath(output_dir)}`.", list(downloaded_accs.intersection(set(all_accessions))), []
        
    total_pending = len(pending_accessions)
    batches = [pending_accessions[i:i + batch_size] for i in range(0, total_pending, batch_size)]
    total_batches = len(batches)
    
    successful_downloads = len(downloaded_accs)
    failed_batches = 0
    
    for batch_idx, batch in enumerate(batches):
        batch_num = batch_idx + 1
        if status_text:
            status_text.info(f"Downloading batch {batch_num}/{total_batches} ({len(batch)} genomes)...")
            
        acc_file = os.path.join(output_dir, f"accessions_batch_{batch_num}.txt")
        with open(acc_file, "w") as f:
            f.write("\n".join(batch))
            
        current_zip = os.path.join(output_dir, f"dataset_batch_{batch_num}.zip")
        cmd = ["./datasets", "download", "genome", "accession", "--inputfile", acc_file, "--filename", current_zip, "--include", "genome"]
        
        try:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                if status_text and line.strip():
                    if "Downloading" in line or "Collecting" in line or "Error" in line:
                        status_text.text(f"Batch {batch_num} Status: {line.strip()}")
            process.wait()
            
            if process.returncode == 0 and os.path.exists(current_zip):
                if status_text:
                    status_text.info(f"Batch {batch_num} downloaded. Extracting files...")
                with zipfile.ZipFile(current_zip, 'r') as zip_ref:
                    zip_ref.extractall(output_dir)
                os.remove(current_zip)
                
                data_dir = os.path.join(output_dir, 'ncbi_dataset', 'data')
                batch_fasta_count = 0
                if os.path.exists(data_dir):
                    for root, dirs, files in os.walk(data_dir):
                        for file in files:
                            if file.endswith('.fna'):
                                src = os.path.join(root, file)
                                dst = os.path.join(output_dir, file)
                                import shutil
                                if not os.path.exists(dst):
                                    shutil.move(src, dst)
                                    batch_fasta_count += 1
                    import shutil
                    shutil.rmtree(os.path.join(output_dir, 'ncbi_dataset'), ignore_errors=True)
                    if os.path.exists(os.path.join(output_dir, 'README.md')):
                        os.remove(os.path.join(output_dir, 'README.md'))
                        
                successful_downloads += batch_fasta_count
            else:
                failed_batches += 1
                if status_text:
                    status_text.warning(f"Batch {batch_num} failed or returned no data (Return code: {process.returncode}). Moving to next batch...")
        except Exception as e:
            failed_batches += 1
            if status_text:
                status_text.warning(f"Exception in batch {batch_num}: {str(e)}")
        finally:
            if os.path.exists(acc_file):
                os.remove(acc_file)
                
        if progress_bar:
            current_progress = successful_downloads / total_genomes
            progress_bar.progress(min(1.0, current_progress))

    final_downloaded = get_downloaded_accessions(output_dir)
    success_list = list(final_downloaded.intersection(set(all_accessions)))
    failed_list = list(set(all_accessions) - set(success_list))

    if failed_batches == 0:
        return True, f"Successfully processed {successful_downloads} genome sequence(s) in folder: `{os.path.abspath(output_dir)}`", success_list, failed_list
    else:
        return False, f"Completed with errors. {failed_batches} batches failed. Successfully saved {successful_downloads}/{total_genomes} genomes.", success_list, failed_list


import shutil
import hashlib

def get_file_hash(filepath):
    """Calculate MD5 hash of a file to detect exact duplicates."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def clean_and_rename_genomes(input_dirs, output_main_dir):
    """
    Clean and organize genomes from multiple directories.
    Uses MD5 hashing to skip physical duplicates.
    Keeps the original filename, avoiding collisions by appending a counter if needed.
    """
    if os.path.exists(output_main_dir):
        shutil.rmtree(output_main_dir)
    os.makedirs(output_main_dir, exist_ok=True)
    report = []
    
    # First pass: collect all files and their hashes
    hash_to_files = {} # hash -> list of file info
    
    for input_dir in input_dirs:
        if not os.path.exists(input_dir):
            continue
        for filename in os.listdir(input_dir):
            if filename.endswith(('.fna', '.fasta', '.fa')):
                filepath = os.path.join(input_dir, filename)
                f_hash = get_file_hash(filepath)
                if f_hash not in hash_to_files:
                    hash_to_files[f_hash] = []
                hash_to_files[f_hash].append({'path': filepath, 'name': filename})

    used_filenames = set()
    
    for f_hash, file_list in hash_to_files.items():
        # Heuristic: Pick the "best" filename among duplicates
        # Priority: 1. Contains GCA/GCF, 2. Shortest name (often cleaner)
        best_file = file_list[0]
        for f in file_list:
            if re.search(r'GC[AF]_\d+\.\d+', f['name']):
                best_file = f
                break
        
        original_name = best_file['name']
        new_filename = original_name
        counter = 1
        
        # Prevent filename collisions if different genomes somehow have the same name
        while new_filename in used_filenames:
            name_part, ext_part = os.path.splitext(original_name)
            new_filename = f"{name_part}_{counter}{ext_part}"
            counter += 1
            
        used_filenames.add(new_filename)
        new_path = os.path.join(output_main_dir, new_filename)
        
        shutil.copy2(best_file['path'], new_path)
        
        report.append({
            "Selected_Source": best_file['name'],
            "Original_Count": len(file_list),
            "New_Name": new_filename,
            "Other_Names": ", ".join([f['name'] for f in file_list if f['name'] != best_file['name']])
        })
                
    return pd.DataFrame(report)

def extract_accession_from_filename(name):
    """Extract GCA/GCF accession from string."""
    match = re.search(r'(GC[AF]_\d+\.\d+)', name)
    return match.group(1) if match else None

def get_tool_path(tool_name):
    """Ensure a tool is installed via Pixi and return its path."""
    import shutil
    tool_path = shutil.which(tool_name)
    if tool_path: return tool_path
    
    global_tool = os.path.join(os.getcwd(), ".pixi", "bin", tool_name)
    if os.path.exists(global_tool):
        return global_tool
        
    return None

def install_tools_via_pixi(status_container):
    """Installs required tools (FastANI, Mash) via Pixi in the background."""
    pixi_home = os.path.join(os.getcwd(), ".pixi")
    os.environ["PIXI_HOME"] = pixi_home
    pixi_bin = os.path.join(pixi_home, "bin", "pixi")
    
    if not os.path.exists(pixi_bin):
        status_container.info("Installing Pixi package manager...")
        cmd_install_pixi = f"export PIXI_HOME={pixi_home} && curl -fsSL https://pixi.sh/install.sh | bash"
        subprocess.run(cmd_install_pixi, shell=True, check=True, capture_output=True)
        
    status_container.info("Installing FastANI and Mash via Pixi...")
    cmd_install_tools = f"export PIXI_HOME={pixi_home} && {pixi_bin} global install fastani mash -c bioconda -c conda-forge"
    subprocess.run(cmd_install_tools, shell=True, check=True, capture_output=True)
    
    return os.path.join(pixi_home, "bin", "fastANI"), os.path.join(pixi_home, "bin", "mash")

def run_fastani_dereplication(input_dir, output_dir, fastani_path, mash_path=None, use_mash=True, ani_threshold=99.9, af_threshold=60.0, threads=4, status_text=None):
    """
    Run FastANI to identify sequence-level duplicates.
    If use_mash is True, uses MASH to pre-filter highly similar genome pairs to speed up FastANI.
    """
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    genomes = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(('.fna', '.fasta', '.fa'))]
    if len(genomes) < 2:
        for g in genomes:
            shutil.copy2(g, os.path.join(output_dir, os.path.basename(g)))
        return pd.DataFrame()
        
    list_file = os.path.join(input_dir, "genome_list.txt")
    with open(list_file, "w") as f:
        f.write("\n".join(genomes))
        
    out_file = os.path.join(input_dir, "fastani_out.txt")
    
    if use_mash and mash_path:
        if status_text:
            status_text.info(f"Step 2.1: Running MASH pre-filtering on {len(genomes)} genomes...")
            
        msh_file = os.path.join(input_dir, "genomes.msh")
        mash_dist_file = os.path.join(input_dir, "mash_dist.txt")
        
        # 1. Sketch
        cmd_sketch = [mash_path, "sketch", "-l", list_file, "-o", msh_file, "-p", str(threads)]
        subprocess.run(cmd_sketch, check=True, capture_output=True)
        
        # 2. Dist
        # MASH distance = 1 - ANI. So distance 0.1 means ~90% ANI.
        # We set threshold loosely (e.g. max dist 0.1) to catch everything that MIGHT pass the FastANI threshold.
        max_dist = 1.0 - (ani_threshold / 100.0) + 0.05 # Add 5% buffer
        cmd_dist = [mash_path, "dist", "-d", str(max_dist), "-p", str(threads), msh_file, msh_file]
        with open(mash_dist_file, "w") as f:
            subprocess.run(cmd_dist, stdout=f, check=True)
            
        # 3. Parse and generate pairs for FastANI
        # Note: MASH dist output sometimes lacks headers or may have extra fields.
        try:
            df_mash = pd.read_csv(mash_dist_file, sep='\t', header=None)
            # Take only the first 5 columns to avoid issues with inconsistent trailing tabs
            df_mash = df_mash.iloc[:, :5]
            df_mash.columns = ['query', 'ref', 'dist', 'pval', 'hashes']
            df_mash = df_mash[df_mash['query'] != df_mash['ref']]
        except Exception as e:
            st.error(f"Error parsing MASH output: {e}")
            df_mash = pd.DataFrame()
        
        if df_mash.empty:
            # No similar pairs found at all! FastANI is not needed.
            for g in genomes:
                shutil.copy2(g, os.path.join(output_dir, os.path.basename(g)))
            if os.path.exists(list_file): os.remove(list_file)
            if os.path.exists(msh_file): os.remove(msh_file)
            if os.path.exists(mash_dist_file): os.remove(mash_dist_file)
            return pd.DataFrame()
            
        # Build connected components from MASH to run FastANI on each cluster
        parent = {}
        def find(i):
            if parent[i] == i: return i
            parent[i] = find(parent[i])
            return parent[i]
        def union(i, j):
            root_i = find(i)
            root_j = find(j)
            if root_i != root_j:
                parent[root_i] = root_j
                
        for g in genomes:
            parent[g] = g
            
        for _, row in df_mash.iterrows():
            if row['query'] in parent and row['ref'] in parent:
                union(row['query'], row['ref'])
                
        cluster_map = {}
        for g in genomes:
            root = find(g)
            if root not in cluster_map:
                cluster_map[root] = []
            cluster_map[root].append(g)
            
        components = list(cluster_map.values())
        
        if status_text:
            status_text.info(f"Step 2.2: Found {len(components)} primary clusters. Running FastANI precision check...")
            
        df_ani_list = []
        for i, comp in enumerate(components):
            if len(comp) > 1:
                comp_list_file = os.path.join(input_dir, f"comp_{i}_list.txt")
                comp_out_file = os.path.join(input_dir, f"comp_{i}_out.txt")
                with open(comp_list_file, "w") as f:
                    f.write("\n".join(comp))
                
                cmd = [fastani_path, "--ql", comp_list_file, "--rl", comp_list_file, "-o", comp_out_file, "-t", str(threads)]
                subprocess.run(cmd, capture_output=True)
                
                if os.path.exists(comp_out_file):
                    try:
                        df_comp = pd.read_csv(comp_out_file, sep='\t', header=None, names=['query', 'ref', 'ani', 'matches', 'total'])
                        df_ani_list.append(df_comp)
                    except Exception:
                        pass
                    os.remove(comp_out_file)
                if os.path.exists(comp_list_file): os.remove(comp_list_file)
                
        if df_ani_list:
            df_ani_full = pd.concat(df_ani_list, ignore_index=True)
            df_ani_full.to_csv(out_file, sep='\t', header=False, index=False)
            
        if os.path.exists(msh_file): os.remove(msh_file)
        if os.path.exists(mash_dist_file): os.remove(mash_dist_file)
        
    else:
        if status_text:
            status_text.info(f"Running FastANI all-vs-all on {len(genomes)} genomes... This might take a while.")
        cmd = [fastani_path, "--ql", list_file, "--rl", list_file, "-o", out_file, "-t", str(threads)]
        subprocess.run(cmd, check=True, capture_output=True)
    
    # Parse FastANI output
    # Format: query, ref, ANI, matches, total
    clusters = [] # list of sets
    if os.path.exists(out_file):
        df_ani = pd.read_csv(out_file, sep='\t', header=None, names=['query', 'ref', 'ani', 'matches', 'total'])
        df_ani['af'] = (df_ani['matches'] / df_ani['total']) * 100
        # Filter for high similarity and high alignment fraction
        df_dup = df_ani[(df_ani['ani'] >= ani_threshold) & (df_ani['af'] >= af_threshold) & (df_ani['query'] != df_ani['ref'])]
        
        # Build connected components
        parent = {}
        def find(i):
            if parent[i] == i: return i
            parent[i] = find(parent[i])
            return parent[i]
        def union(i, j):
            root_i = find(i)
            root_j = find(j)
            if root_i != root_j:
                parent[root_i] = root_j
                
        for g in genomes:
            parent[g] = g
            
        for _, row in df_dup.iterrows():
            union(row['query'], row['ref'])
            
        cluster_map = {}
        for g in genomes:
            root = find(g)
            if root not in cluster_map:
                cluster_map[root] = []
            cluster_map[root].append(g)
            
        report = []
        for root, members in cluster_map.items():
            # Sort to pick the best representative (e.g. shortest name or containing GCA)
            members.sort(key=lambda x: (0 if re.search(r'GC[AF]_', x) else 1, len(x)))
            rep = members[0]
            shutil.copy2(rep, os.path.join(output_dir, os.path.basename(rep)))
            
            for dup in members[1:]:
                report.append({
                    "Representative": os.path.basename(rep),
                    "Duplicate_Found": os.path.basename(dup),
                    "ANI": df_dup[(df_dup['query']==dup) & (df_dup['ref']==rep)]['ani'].max() if not df_dup[(df_dup['query']==dup) & (df_dup['ref']==rep)].empty else df_dup[(df_dup['query']==rep) & (df_dup['ref']==dup)]['ani'].max()
                })
                
        if os.path.exists(list_file): os.remove(list_file)
        if os.path.exists(out_file): os.remove(out_file)
        
        return pd.DataFrame(report)
    return pd.DataFrame()

def run_fastani_comparison(my_genomes_dir, gtdb_genomes_dir, fastani_path, mash_path=None, use_mash=True, ani_threshold=95.0, af_threshold=60.0, threads=4, status_text=None):
    """
    Compare user dataset vs GTDB downloaded dataset.
    If use_mash is True, uses MASH to pre-filter to only compute FastANI on pairs with similar distance.
    """
    my_genomes = [os.path.join(my_genomes_dir, f) for f in os.listdir(my_genomes_dir) if f.endswith(('.fna', '.fasta', '.fa'))]
    gtdb_genomes = [os.path.join(gtdb_genomes_dir, f) for f in os.listdir(gtdb_genomes_dir) if f.endswith(('.fna', '.fasta', '.fa'))]
    
    if not my_genomes or not gtdb_genomes:
        return None
        
    my_list = os.path.join(my_genomes_dir, "my_query_list.txt")
    with open(my_list, "w") as f: f.write("\n".join(my_genomes))
        
    gtdb_list = os.path.join(gtdb_genomes_dir, "gtdb_ref_list.txt")
    with open(gtdb_list, "w") as f: f.write("\n".join(gtdb_genomes))
        
    out_file = os.path.join(my_genomes_dir, "fastani_vs_gtdb.out")
    
    if use_mash and mash_path:
        if status_text:
            status_text.info(f"Running MASH pre-filtering: {len(my_genomes)} My Genomes VS {len(gtdb_genomes)} GTDB Genomes...")
            
        my_msh = os.path.join(my_genomes_dir, "my.msh")
        gtdb_msh = os.path.join(gtdb_genomes_dir, "gtdb.msh")
        mash_dist_file = os.path.join(my_genomes_dir, "mash_vs_gtdb.txt")
        
        # Sketch both
        subprocess.run([mash_path, "sketch", "-l", my_list, "-o", my_msh, "-p", str(threads)], check=True, capture_output=True)
        subprocess.run([mash_path, "sketch", "-l", gtdb_list, "-o", gtdb_msh, "-p", str(threads)], check=True, capture_output=True)
        
        # Dist
        max_dist = 1.0 - (ani_threshold / 100.0) + 0.05
        cmd_dist = [mash_path, "dist", "-d", str(max_dist), "-p", str(threads), gtdb_msh, my_msh]
        with open(mash_dist_file, "w") as f:
            subprocess.run(cmd_dist, stdout=f, check=True)
            
        # Parse
        try:
            df_mash = pd.read_csv(mash_dist_file, sep='\t', header=None)
            df_mash = df_mash.iloc[:, :5]
            df_mash.columns = ['ref', 'query', 'dist', 'pval', 'hashes']
        except Exception as e:
            st.error(f"Error parsing MASH comparison output: {e}")
            df_mash = pd.DataFrame()
        
        if df_mash.empty:
            if os.path.exists(my_list): os.remove(my_list)
            if os.path.exists(gtdb_list): os.remove(gtdb_list)
            if os.path.exists(my_msh): os.remove(my_msh)
            if os.path.exists(gtdb_msh): os.remove(gtdb_msh)
            if os.path.exists(mash_dist_file): os.remove(mash_dist_file)
            return pd.DataFrame()
            
        unique_queries = df_mash['query'].unique()
        unique_refs = df_mash['ref'].unique()
        
        filtered_my_list = os.path.join(my_genomes_dir, "my_query_filtered_list.txt")
        filtered_gtdb_list = os.path.join(gtdb_genomes_dir, "gtdb_ref_filtered_list.txt")
        
        with open(filtered_my_list, "w") as f:
            f.write("\n".join(unique_queries))
        with open(filtered_gtdb_list, "w") as f:
            f.write("\n".join(unique_refs))
            
        if status_text:
            status_text.info(f"Running FastANI on {len(unique_queries)} vs {len(unique_refs)} pre-filtered genomes...")
        
        cmd = [fastani_path, "--ql", filtered_my_list, "--rl", filtered_gtdb_list, "-o", out_file, "-t", str(threads)]
        subprocess.run(cmd, capture_output=True)
        
        if os.path.exists(my_msh): os.remove(my_msh)
        if os.path.exists(gtdb_msh): os.remove(gtdb_msh)
        if os.path.exists(mash_dist_file): os.remove(mash_dist_file)
        if os.path.exists(filtered_my_list): os.remove(filtered_my_list)
        if os.path.exists(filtered_gtdb_list): os.remove(filtered_gtdb_list)
        
    else:
        if status_text:
            status_text.info(f"Running FastANI comparison: {len(my_genomes)} My Genomes VS {len(gtdb_genomes)} GTDB Genomes...")
        cmd = [fastani_path, "--ql", my_list, "--rl", gtdb_list, "-o", out_file, "-t", str(threads)]
        subprocess.run(cmd, check=True, capture_output=True)
    
    if os.path.exists(out_file):
        try:
            df_ani = pd.read_csv(out_file, sep='\t', header=None, names=['query', 'ref', 'ani', 'matches', 'total'])
            df_ani['query_name'] = df_ani['query'].apply(os.path.basename)
            df_ani['ref_name'] = df_ani['ref'].apply(os.path.basename)
            df_ani['af'] = (df_ani['matches'] / df_ani['total']) * 100
            
            # Filter by threshold
            df_match = df_ani[(df_ani['ani'] >= ani_threshold) & (df_ani['af'] >= af_threshold)]
            
            if os.path.exists(my_list): os.remove(my_list)
            if os.path.exists(gtdb_list): os.remove(gtdb_list)
            
            return df_match
        except Exception:
            return pd.DataFrame()
            
    return pd.DataFrame()

st.title("GTDB & NCBI Genomes Toolkit")

df = load_data()
if df.empty:
    st.error("No data found in 'gtdb_data' directory. Please run the download script first.")
    st.stop()

versions = sorted(df['Version'].unique())

tab1, tab2, tab3, tab4 = st.tabs(["🔍 Single Version Explorer", "📊 Version Comparison", "📥 Custom Download", "📦 Dataset Updater"])

with tab1:
    st.header("Single Version Explorer")
    st.markdown("Explore taxonomy and genome information for a specific GTDB release.")
    sel_version = st.selectbox("Select GTDB Version", versions, index=len(versions)-1)
    
    df_single = df[df['Version'] == sel_version]
    st.write(f"**Total Genomes in R{sel_version}:** {len(df_single)}")
    
    search_single = st.text_input("Enter a taxonomic group to search (e.g., c__Bathyarchaeia):", "c__Bathyarchaeia", key="search_single")
    if search_single:
        df_res = df_single[df_single['Taxonomy'].str.contains(search_single, na=False)].copy()
        rep_count = df_res['Species'].nunique()
        st.markdown(f"**Found {len(df_res)} total genomes** and **{rep_count} representative genomes (unique species)** matching `{search_single}`.")
        
        # Mark one representative per species
        df_res['Is_Representative'] = 'No'
        rep_idx = df_res.drop_duplicates(subset=['Species']).index
        df_res.loc[rep_idx, 'Is_Representative'] = 'Yes'
        
        st.info("Note: Since we only load taxonomy data, one genome per unique species is automatically marked as the representative.")
        
        def highlight_reps(row):
            return ['background-color: rgba(40, 167, 69, 0.2)'] * len(row) if row['Is_Representative'] == 'Yes' else [''] * len(row)
            
        st.dataframe(df_res[['Genome_ID', 'Is_Representative', 'Taxonomy']].style.apply(highlight_reps, axis=1))
        
        col_dl1, col_dl2 = st.columns(2)
        with col_dl1:
            if st.button("Download ALL genomes via NCBI Datasets", key="btn_single_dl_all"):
                if ensure_datasets_cli():
                    st_cont = st.empty()
                    p_bar = st.progress(0)
                    safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_single)
                    out_dir = f"ncbi_downloads/{safe_name}_R{sel_version}_all"
                    succ, msg, s_list, f_list = download_genomes(df_res['Genome_ID'].tolist(), out_dir, progress_bar=p_bar, status_text=st_cont)
                    if succ:
                        st_cont.success(msg)
                    else:
                        st_cont.warning(msg)
                else:
                    st.error("NCBI Datasets CLI could not be installed.")
        with col_dl2:
            if st.button("Download ONLY Representatives via NCBI Datasets", key="btn_single_dl_reps"):
                if ensure_datasets_cli():
                    st_cont = st.empty()
                    p_bar = st.progress(0)
                    safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_single)
                    out_dir = f"ncbi_downloads/{safe_name}_R{sel_version}_reps"
                    reps_list = df_res[df_res['Is_Representative'] == 'Yes']['Genome_ID'].tolist()
                    succ, msg, s_list, f_list = download_genomes(reps_list, out_dir, progress_bar=p_bar, status_text=st_cont)
                    if succ:
                        st_cont.success(msg)
                    else:
                        st_cont.warning(msg)
                else:
                    st.error("NCBI Datasets CLI could not be installed.")

with tab2:
    st.header("Version Comparison")
    st.markdown("Compare taxonomy changes between two GTDB releases.")
    selected_versions = st.multiselect("Select Versions to Compare", versions, default=versions[-2:] if len(versions)>=2 else versions)
    selected_versions = sorted(selected_versions)
    
    if len(selected_versions) < 2:
        st.warning("Please select at least 2 versions to compare.")
    else:
        search_term = st.text_input("Enter a taxonomic group to analyze:", "c__Bathyarchaeia", key="search_comp")
        if search_term:
            v1 = selected_versions[-2]
            v2 = selected_versions[-1]
            st.markdown(f"**Comparing R{v1} vs R{v2}**")
            
            df_v1 = df[(df['Version'] == v1) & (df['Taxonomy'].str.contains(search_term, na=False))]
            df_v2 = df[(df['Version'] == v2) & (df['Taxonomy'].str.contains(search_term, na=False))]
            
            st.markdown("##### 🧬 Total Genomes")
            col1, col2, col3 = st.columns(3)
            col1.metric("Genomes in R" + str(v1), len(df_v1))
            col2.metric("Genomes in R" + str(v2), len(df_v2))
            col3.metric("Net Change", len(df_v2) - len(df_v1))
            
            st.markdown("##### 👑 Representative Genomes (Unique Species)")
            rep_v1 = df_v1['Species'].nunique()
            rep_v2 = df_v2['Species'].nunique()
            col4, col5, col6 = st.columns(3)
            col4.metric("Reps in R" + str(v1), rep_v1)
            col5.metric("Reps in R" + str(v2), rep_v2)
            col6.metric("Rep Net Change", rep_v2 - rep_v1)
            
            genomes_v1 = set(df_v1['Genome_ID'])
            genomes_v2 = set(df_v2['Genome_ID'])
            new_genomes = genomes_v2 - genomes_v1
            
            st.write(f"- **Newly added genomes:** {len(new_genomes)}")
            st.write(f"- **Genomes no longer in this group:** {len(genomes_v1 - genomes_v2)}")
            
            # --- New Detailed Comparison Section ---
            st.markdown("### Detailed Comparison Details")
            
            diff_tab1, diff_tab2, diff_tab3 = st.tabs(["🆕 Added Genomes", "❌ Removed/Moved Genomes", "🔄 Taxonomy Updates"])
            
            with diff_tab1:
                if new_genomes:
                    st.write(f"Showing {len(new_genomes)} genomes added in R{v2}:")
                    st.dataframe(df_v2[df_v2['Genome_ID'].isin(new_genomes)][['Genome_ID', 'Taxonomy']])
                else:
                    st.info("No new genomes added in this version.")
                    
            with diff_tab2:
                removed_genomes = genomes_v1 - genomes_v2
                if removed_genomes:
                    st.write(f"Showing {len(removed_genomes)} genomes removed from this group in R{v2}:")
                    # Find where they went in R2
                    lost_df = df_v1[df_v1['Genome_ID'].isin(removed_genomes)][['Genome_ID', 'Taxonomy']]
                    where_went = df[(df['Version'] == v2) & (df['Genome_ID'].isin(removed_genomes))][['Genome_ID', 'Taxonomy']]
                    where_went = where_went.rename(columns={'Taxonomy': f'New_Taxonomy_R{v2}'})
                    lost_compare = lost_df.merge(where_went, on='Genome_ID', how='left')
                    st.dataframe(lost_compare)
                else:
                    st.info("No genomes were removed or moved.")
                    
            with diff_tab3:
                retained_genomes = genomes_v1 & genomes_v2
                df_retained_v1 = df_v1[df_v1['Genome_ID'].isin(retained_genomes)][['Genome_ID', 'Taxonomy']].set_index('Genome_ID')
                df_retained_v2 = df_v2[df_v2['Genome_ID'].isin(retained_genomes)][['Genome_ID', 'Taxonomy']].set_index('Genome_ID')
                retained_compare = df_retained_v1.join(df_retained_v2, lsuffix=f'_R{v1}', rsuffix=f'_R{v2}')
                changed_tax = retained_compare[retained_compare[f'Taxonomy_R{v1}'] != retained_compare[f'Taxonomy_R{v2}']]
                
                if not changed_tax.empty:
                    st.write(f"Showing {len(changed_tax)} genomes with updated taxonomy strings in R{v2}:")
                    st.dataframe(changed_tax)
                else:
                    st.info("No taxonomy updates found for retained genomes.")
            # --- End of Detailed Comparison Section ---

            dl_col1, dl_col2, dl_col3, dl_col4 = st.columns(4)
            with dl_col1:
                if st.button(f"Download ALL Genomes ({len(df_v2)}) in R{v2}", key="dl_all"):
                    if ensure_datasets_cli():
                        st_cont = st.empty()
                        p_bar = st.progress(0)
                        safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_term)
                        out_dir = f"ncbi_downloads/{safe_name}_R{v2}_all"
                        succ, msg, s_list, f_list = download_genomes(df_v2['Genome_ID'].tolist(), out_dir, progress_bar=p_bar, status_text=st_cont)
                        if succ:
                            st_cont.success(msg)
                        else:
                            st_cont.warning(msg)
                    else:
                        st.error("NCBI Datasets CLI could not be installed.")
            with dl_col2:
                if len(new_genomes) > 0 and st.button(f"Download NEW Genomes ({len(new_genomes)}) in R{v2}", key="dl_new"):
                    if ensure_datasets_cli():
                        st_cont = st.empty()
                        p_bar = st.progress(0)
                        safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_term)
                        out_dir = f"ncbi_downloads/{safe_name}_R{v2}_new"
                        succ, msg, s_list, f_list = download_genomes(list(new_genomes), out_dir, progress_bar=p_bar, status_text=st_cont)
                        if succ:
                            st_cont.success(msg)
                        else:
                            st_cont.warning(msg)
                    else:
                        st.error("NCBI Datasets CLI could not be installed.")
            with dl_col3:
                reps_v2 = df_v2.drop_duplicates(subset=['Species'])['Genome_ID'].tolist()
                if st.button(f"Download ALL Reps ({len(reps_v2)}) in R{v2}", key="dl_all_reps"):
                    if ensure_datasets_cli():
                        st_cont = st.empty()
                        p_bar = st.progress(0)
                        safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_term)
                        out_dir = f"ncbi_downloads/{safe_name}_R{v2}_all_reps"
                        succ, msg, s_list, f_list = download_genomes(reps_v2, out_dir, progress_bar=p_bar, status_text=st_cont)
                        if succ:
                            st_cont.success(msg)
                        else:
                            st_cont.warning(msg)
                    else:
                        st.error("NCBI Datasets CLI could not be installed.")
            with dl_col4:
                new_species = set(df_v2['Species']) - set(df_v1['Species'])
                new_reps = df_v2[df_v2['Species'].isin(new_species)].drop_duplicates(subset=['Species'])['Genome_ID'].tolist()
                if len(new_reps) > 0 and st.button(f"Download NEW Reps ({len(new_reps)}) in R{v2}", key="dl_new_reps"):
                    if ensure_datasets_cli():
                        st_cont = st.empty()
                        p_bar = st.progress(0)
                        safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_term)
                        out_dir = f"ncbi_downloads/{safe_name}_R{v2}_new_reps"
                        succ, msg, s_list, f_list = download_genomes(new_reps, out_dir, progress_bar=p_bar, status_text=st_cont)
                        if succ:
                            st_cont.success(msg)
                        else:
                            st_cont.warning(msg)
                    else:
                        st.error("NCBI Datasets CLI could not be installed.")

with tab3:
    st.header("Custom Download")
    st.markdown("Enter custom Genome IDs to download them directly via NCBI Datasets.")
    
    custom_ids_text = st.text_area("Enter Genome IDs (one per line, e.g., GCF_000000000.1, RS_GCF_000979855.1):")
    
    if st.button("Start Download"):
        ids = [i.strip() for i in custom_ids_text.split('\n') if i.strip()]
        if not ids:
            st.error("Please enter at least one Genome ID.")
        else:
            if ensure_datasets_cli():
                st_cont = st.empty()
                p_bar = st.progress(0)
                out_dir = "ncbi_downloads/custom_download"
                
                st_cont.info("Starting download process...")
                succ, msg, s_list, f_list = download_genomes(ids, out_dir, progress_bar=p_bar, status_text=st_cont)
                
                if succ:
                    st_cont.success(msg)
                else:
                    st_cont.warning(msg)
                    
                # Generate Report
                report_file = "ncbi_downloads/custom_download_report.csv"
                report_df = pd.DataFrame({
                    "Requested_ID": ids,
                    "Normalized_Accession": [g[3:] if g.startswith(("RS_", "GB_")) else g for g in ids],
                })
                report_df["Status"] = report_df["Normalized_Accession"].apply(lambda x: "Success" if x in s_list else "Failed")
                report_df.to_csv(report_file, index=False)
                
                st.success(f"🎉 Download finished! Detailed report saved to `{report_file}`.")
            else:
                st.error("NCBI Datasets CLI could not be installed.")

with tab4:
    st.header("Local Dataset vs GTDB Species-Level Pipeline")
    st.markdown("""
    **What are `matches` and `total`?**
    In FastANI, genomes are conceptually divided into sequence fragments (usually 3kb long). 
    - `matches`: The number of fragments from the query genome that successfully aligned with the reference genome.
    - `total`: The total number of sequence fragments in the query genome.
    - **Alignment Fraction (AF)**: Calculated as `(matches / total) * 100`. It represents the coverage of the alignment.
    
    **Pipeline Workflow:**
    This one-click pipeline performs a strict species-level comparison:
    1. **Physical Cleaning**: Deduplicates identical local files based on MD5 hashing.
    2. **Local Dereplication**: Clusters your local genomes into **Unique Species** (Local Representatives) using FastANI & MASH.
    3. **GTDB Auto-Fetch**: Automatically identifies and downloads the **Species Representatives** for your specified GTDB Taxon.
    4. **Species vs Species Comparison**: Compares your Local Species against the GTDB Species to reveal overlapping and novel species.
    """)
    
    with st.expander("Pipeline Configuration", expanded=True):
        input_dirs_raw = st.text_area("1. Local Genome Folders (one path per line):", 
                                     help="Example: /Users/user/genomes/bathy_v1\n/Users/user/genomes/bathy_v2")
        output_name = st.text_input("2. Dataset Name (for saving results):", "Bathyarchaeia_Species_Comparison")
        
        st.markdown("##### 🧬 Local Species Definition (Dereplication)")
        col_d1, col_d2 = st.columns(2)
        with col_d1:
            derep_ani_thresh = st.slider("Local Species ANI Threshold (%)", min_value=80.0, max_value=100.0, value=95.0, step=0.1, help="Typically 95% for species-level clustering.")
        with col_d2:
            derep_af_thresh = st.slider("Local Species AF Threshold (%)", min_value=10.0, max_value=100.0, value=65.0, step=1.0, help="Typically 65% coverage for robust species definition.")
            
        st.markdown("##### 🌍 GTDB Target (For Comparison)")
        col_g1, col_g2 = st.columns(2)
        with col_g1:
            gtdb_comp_version = st.selectbox("GTDB Version", versions, index=len(versions)-1, key="gtdb_comp_v")
        with col_g2:
            gtdb_comp_taxon = st.text_input("GTDB Taxon to compare against:", "c__Bathyarchaeia", key="gtdb_comp_t")
            
        st.markdown("##### ⚖️ Comparison Matching Rules")
        col_c1, col_c2 = st.columns(2)
        with col_c1:
            comp_ani_thresh = st.slider("Comparison Match ANI Threshold (%)", min_value=80.0, max_value=100.0, value=95.0, step=0.1)
        with col_c2:
            comp_af_thresh = st.slider("Comparison Match AF Threshold (%)", min_value=10.0, max_value=100.0, value=65.0, step=1.0)
            
        st.markdown("##### ⚙️ Compute Settings")
        col_m1, col_m2 = st.columns(2)
        with col_m1:
            threads = st.number_input("Threads to use for FastANI & MASH", min_value=1, max_value=64, value=4)
        with col_m2:
            st.markdown("<br>", unsafe_allow_html=True)
            use_mash_pipe = st.checkbox("Use MASH pre-filtering (Highly Recommended)", value=True)
            
    if st.button("Run Full Species-Level Analysis"):
        input_dirs = [d.strip() for d in input_dirs_raw.split('\n') if d.strip()]
        if not input_dirs:
            st.error("Please provide at least one local input directory.")
        elif not gtdb_comp_taxon:
            st.error("Please specify a GTDB Taxon.")
        else:
            status_container = st.empty()
            
            # --- Sub-step A: Fetch GTDB Representatives ---
            status_container.info(f"Step A: Retrieving GTDB Representatives for {gtdb_comp_taxon} (R{gtdb_comp_version})...")
            df_gtdb = df[(df['Version'] == gtdb_comp_version) & (df['Taxonomy'].str.contains(gtdb_comp_taxon, na=False))]
            if df_gtdb.empty:
                st.error(f"No genomes found for {gtdb_comp_taxon} in GTDB R{gtdb_comp_version}.")
                st.stop()
                
            df_gtdb['Is_Representative'] = 'No'
            rep_idx = df_gtdb.drop_duplicates(subset=['Species']).index
            df_gtdb.loc[rep_idx, 'Is_Representative'] = 'Yes'
            reps_list = df_gtdb[df_gtdb['Is_Representative'] == 'Yes']['Genome_ID'].tolist()
            
            safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', gtdb_comp_taxon)
            gtdb_dir = f"ncbi_downloads/{safe_name}_R{gtdb_comp_version}_reps"
            os.makedirs(gtdb_dir, exist_ok=True)
            
            downloaded = get_downloaded_accessions(gtdb_dir)
            missing = [acc for acc in reps_list if acc.replace("RS_", "").replace("GB_", "") not in downloaded and acc not in downloaded]
            
            if missing:
                status_container.info(f"Step A: Downloading {len(missing)} missing GTDB Representative Genomes via NCBI Datasets...")
                if ensure_datasets_cli():
                    p_bar = st.progress(0)
                    succ, msg, s_list, f_list = download_genomes(reps_list, gtdb_dir, progress_bar=p_bar, status_text=status_container)
                    p_bar.empty()
                    if not succ:
                        st.warning(msg)
                else:
                    st.error("NCBI Datasets CLI could not be installed.")
                    st.stop()
            else:
                status_container.info(f"Step A: All {len(reps_list)} GTDB Representative Genomes are already downloaded locally.")
                
            # --- Sub-step B: Physical Cleaning ---
            combined_out = f"local_datasets/{output_name}/1_raw_combined"
            derep_out = f"local_datasets/{output_name}/2_dereplicated_genomes"
            status_container.info("Step B: Physical Hashing and Cleaning of local genomes...")
            report_df = clean_and_rename_genomes(input_dirs, combined_out)
            
            if report_df.empty:
                st.warning("No genome files (.fna, .fasta, .fa) found in the provided directories.")
                st.stop()
                
            # --- Sub-step C: Local Dereplication ---
            fastani_path = get_tool_path("fastANI")
            mash_path = get_tool_path("mash")
            
            if not fastani_path or not mash_path:
                with st.spinner("Installing Pixi, FastANI, and Mash..."):
                    try:
                        fastani_path, mash_path = install_tools_via_pixi(status_container)
                    except Exception as e:
                        st.error(f"Failed to install tools: {e}")
                        st.stop()
                        
            status_container.info(f"Step C: Dereplicating local genomes into Species (ANI {derep_ani_thresh}%, AF {derep_af_thresh}%)...")
            dup_report = run_fastani_dereplication(
                combined_out, derep_out, fastani_path, mash_path=mash_path, use_mash=use_mash_pipe, 
                ani_threshold=derep_ani_thresh, af_threshold=derep_af_thresh, threads=threads, status_text=status_container
            )
            
            # --- Sub-step D: Compare Local Species vs GTDB Species ---
            status_container.info("Step D: Running Species vs Species Comparison...")
            df_match = run_fastani_comparison(
                derep_out, gtdb_dir, fastani_path, mash_path=mash_path, use_mash=use_mash_pipe, 
                ani_threshold=comp_ani_thresh, af_threshold=comp_af_thresh, threads=threads, status_text=status_container
            )
            
            status_container.empty()
            
            # --- Sub-step E: Render Results ---
            if df_match is not None:
                st.success("🎉 Species-Level Analysis Complete!")
                
                local_reps = set([f for f in os.listdir(derep_out) if f.endswith(('.fna', '.fasta', '.fa'))])
                gtdb_reps = set([f for f in os.listdir(gtdb_dir) if f.endswith(('.fna', '.fasta', '.fa'))])
                
                matched_local = set(df_match['query_name']) if not df_match.empty else set()
                matched_gtdb = set(df_match['ref_name']) if not df_match.empty else set()
                
                novel_local = local_reps - matched_local
                missing_gtdb = gtdb_reps - matched_gtdb
                
                st.markdown("### Species Comparison Summary")
                st.markdown("This strictly compares **Unique Local Species** against **GTDB Representative Species**.")
                
                col_r1, col_r2, col_r3, col_r4 = st.columns(4)
                col_r1.metric("Local Unique Species", len(local_reps))
                col_r2.metric("GTDB Unique Species", len(gtdb_reps))
                col_r3.metric("Shared Species (Overlap)", len(matched_local))
                col_r4.metric("Novel Local Species", len(novel_local))
                
                c_tab1, c_tab2, c_tab3, c_tab4 = st.tabs(["🏠 Novel Local Species", "🌍 Missing GTDB Species", "🤝 Shared Species", "🗑️ Local Duplicates Removed"])
                
                with c_tab1:
                    st.write(f"**{len(novel_local)}** novel species found in your dataset that are not in GTDB R{gtdb_comp_version}:")
                    st.dataframe(pd.DataFrame({"Novel_Local_Species_File": list(novel_local)}))
                    
                with c_tab2:
                    st.write(f"**{len(missing_gtdb)}** species in GTDB R{gtdb_comp_version} that you haven't collected locally:")
                    st.dataframe(pd.DataFrame({"Missing_GTDB_Species_File": list(missing_gtdb)}))
                    
                with c_tab3:
                    st.write(f"**{len(df_match)}** match relationships between your local species and GTDB species:")
                    if not df_match.empty:
                        st.dataframe(df_match[['query_name', 'ref_name', 'ani', 'af', 'matches', 'total']])
                    else:
                        st.info("No overlap found.")
                        
                with c_tab4:
                    st.write("These genomes were identified as duplicates and merged into their respective Local Species Representatives:")
                    total_files = report_df['Original_Count'].sum()
                    st.write(f"- Total Initial Files: {total_files}")
                    st.write(f"- Exact MD5 Duplicates Removed: {total_files - len(report_df)}")
                    st.write(f"- Sequence Duplicates Removed (ANI/AF): {len(dup_report) if not dup_report.empty else 0}")
                    if not dup_report.empty:
                        st.dataframe(dup_report)
            else:
                st.error("Comparison failed or one of the directories contained no fasta files.")

