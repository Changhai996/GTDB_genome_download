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

def send_email_report(sender_email, sender_pwd, smtp_server, smtp_port, receiver_email, report_file, summary_text):
    msg = MIMEMultipart()
    msg['Subject'] = 'NCBI Datasets Download Report'
    msg['From'] = sender_email
    msg['To'] = receiver_email
    
    msg.attach(MIMEText(summary_text, 'plain'))
    
    if os.path.exists(report_file):
        with open(report_file, 'rb') as f:
            part = MIMEApplication(f.read(), Name=os.path.basename(report_file))
            part['Content-Disposition'] = f'attachment; filename="{os.path.basename(report_file)}"'
            msg.attach(part)
        
    try:
        if smtp_port == 465:
            with smtplib.SMTP_SSL(smtp_server, smtp_port) as server:
                server.login(sender_email, sender_pwd)
                server.send_message(msg)
        else:
            with smtplib.SMTP(smtp_server, smtp_port) as server:
                server.starttls()
                server.login(sender_email, sender_pwd)
                server.send_message(msg)
        return True, "Email sent successfully."
    except Exception as e:
        return False, f"Failed to send email: {e}"

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

def get_fastani_path():
    """Ensure FastANI is installed via Pixi and return its path."""
    import shutil
    fastani = shutil.which("fastANI")
    if fastani: return fastani
    
    # Check project local pixi installation
    local_pixi = os.path.join(os.getcwd(), ".pixi", "bin", "pixi")
    local_fastani = os.path.join(os.getcwd(), ".pixi", "envs", "default", "bin", "fastANI")
    global_fastani = os.path.join(os.getcwd(), ".pixi", "bin", "fastANI")
    
    if os.path.exists(global_fastani):
        return global_fastani
        
    return None

def install_fastani(status_container):
    """Installs FastANI via Pixi in the background."""
    pixi_home = os.path.join(os.getcwd(), ".pixi")
    os.environ["PIXI_HOME"] = pixi_home
    pixi_bin = os.path.join(pixi_home, "bin", "pixi")
    
    if not os.path.exists(pixi_bin):
        status_container.info("Installing Pixi package manager...")
        cmd_install_pixi = f"export PIXI_HOME={pixi_home} && curl -fsSL https://pixi.sh/install.sh | bash"
        subprocess.run(cmd_install_pixi, shell=True, check=True, capture_output=True)
        
    status_container.info("Installing FastANI via Pixi...")
    cmd_install_fastani = f"export PIXI_HOME={pixi_home} && {pixi_bin} global install fastani -c bioconda -c conda-forge"
    subprocess.run(cmd_install_fastani, shell=True, check=True, capture_output=True)
    
    return os.path.join(pixi_home, "bin", "fastANI")

def run_fastani_dereplication(input_dir, output_dir, fastani_path, ani_threshold=99.9, af_threshold=60.0, threads=4, status_text=None):
    """
    Run FastANI all-vs-all to identify sequence-level duplicates.
    Keep one representative per cluster and move to output_dir.
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
        
    out_file = os.path.join(input_dir, "fastani_all_vs_all.out")
    
    if status_text:
        status_text.info(f"Running FastANI all-vs-all on {len(genomes)} genomes... This might take a few minutes.")
        
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

def run_fastani_comparison(my_genomes_dir, gtdb_genomes_dir, fastani_path, ani_threshold=95.0, af_threshold=60.0, threads=4, status_text=None):
    """
    Compare user dataset vs GTDB downloaded dataset.
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

tab1, tab2, tab3, tab4 = st.tabs(["🔍 Single Version Explorer", "📊 Version Comparison", "📥 Custom Download & Email", "📦 Dataset Updater"])

with tab1:
    st.header("Single Version Explorer")
    st.markdown("Explore taxonomy and genome information for a specific GTDB release.")
    sel_version = st.selectbox("Select GTDB Version", versions, index=len(versions)-1)
    
    df_single = df[df['Version'] == sel_version]
    st.write(f"**Total Genomes in R{sel_version}:** {len(df_single)}")
    
    search_single = st.text_input("Enter a taxonomic group to search (e.g., c__Bathyarchaeia):", "c__Bathyarchaeia", key="search_single")
    if search_single:
        df_res = df_single[df_single['Taxonomy'].str.contains(search_single, na=False)]
        rep_count = df_res['Species'].nunique()
        st.markdown(f"**Found {len(df_res)} total genomes** and **{rep_count} representative genomes (unique species)** matching `{search_single}`.")
        st.dataframe(df_res[['Genome_ID', 'Taxonomy']])
        
        if st.button("Download these genomes via NCBI Datasets", key="btn_single_dl"):
            if ensure_datasets_cli():
                st_cont = st.empty()
                p_bar = st.progress(0)
                safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', search_single)
                out_dir = f"ncbi_downloads/{safe_name}_R{sel_version}"
                succ, msg, s_list, f_list = download_genomes(df_res['Genome_ID'].tolist(), out_dir, progress_bar=p_bar, status_text=st_cont)
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

            dl_col1, dl_col2 = st.columns(2)
            with dl_col1:
                if st.button(f"Download ALL ({len(df_v2)}) in R{v2}", key="dl_all"):
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
            with dl_col2:
                if len(new_genomes) > 0 and st.button(f"Download NEW ({len(new_genomes)}) in R{v2}", key="dl_new"):
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

with tab3:
    st.header("Custom Download & Email Report")
    st.markdown("Enter custom Genome IDs to download, and automatically send a report to your email.")
    
    custom_ids_text = st.text_area("Enter Genome IDs (one per line, e.g., GCF_000000000.1, RS_GCF_000979855.1):")
    
    st.subheader("Email Configuration")
    st.info("To send an email, please provide a valid sender email and its SMTP password/authorization code.")
    
    col_em1, col_em2 = st.columns(2)
    with col_em1:
        sender_email = st.text_input("Sender Email (e.g., your_email@163.com)")
        sender_pwd = st.text_input("Sender Password / Auth Code", type="password")
    with col_em2:
        smtp_server = st.text_input("SMTP Server (e.g., smtp.163.com)", "smtp.163.com")
        smtp_port = st.number_input("SMTP Port (usually 465 for SSL or 25/587)", value=465)
        
    receiver_email = st.text_input("Receiver Email", "changhaiduan@163.com")
    
    if st.button("Start Download & Send Email"):
        ids = [i.strip() for i in custom_ids_text.split('\n') if i.strip()]
        if not ids:
            st.error("Please enter at least one Genome ID.")
        elif not (sender_email and sender_pwd and smtp_server and receiver_email):
            st.error("Please fill in all Email Configuration fields to receive the report.")
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
                
                # Prepare Email Summary
                summary = f"NCBI Datasets Download Summary:\n\n"
                summary += f"Total Requested: {len(ids)}\n"
                summary += f"Successfully Downloaded: {len(s_list)}\n"
                summary += f"Failed to Download: {len(f_list)}\n\n"
                summary += "Please find the detailed status report attached."
                
                st_cont.info("Sending email report...")
                email_succ, email_msg = send_email_report(sender_email, sender_pwd, smtp_server, smtp_port, receiver_email, report_file, summary)
                
                if email_succ:
                    st.success("🎉 Download finished and email sent successfully!")
                else:
                    st.error(f"Download finished, but failed to send email: {email_msg}")
            else:
                st.error("NCBI Datasets CLI could not be installed.")

with tab4:
    st.header("Local Dataset Updater & FastANI Dereplication")
    st.markdown("""
    Integrate genomes from multiple local folders, automatically run **FastANI** to identify and remove sequence-level duplicates, 
    and compare your unique dataset directly against downloaded GTDB sequences.
    """)
    
    with st.expander("Step 1: Configure Input Directories", expanded=True):
        input_dirs_raw = st.text_area("Enter paths to local genome folders (one per line):", 
                                     help="Example: /Users/user/genomes/bathy_v1\n/Users/user/genomes/bathy_v2")
        output_name = st.text_input("New Dataset Name:", "Bathyarchaeia_Combined")
        derep_ani_thresh = st.slider("Dereplication ANI Threshold (%)", min_value=99.0, max_value=100.0, value=99.9, step=0.1, help="Genomes with ANI >= this threshold will be considered identical duplicates.")
        derep_af_thresh = st.slider("Dereplication AF Threshold (%)", min_value=10.0, max_value=100.0, value=60.0, step=1.0, help="Genomes must also have an Alignment Fraction (AF) >= this threshold to be considered duplicates.")
        threads = st.number_input("Threads to use for FastANI", min_value=1, max_value=64, value=4)
        
    if st.button("Integrate & Dereplicate"):
        input_dirs = [d.strip() for d in input_dirs_raw.split('\n') if d.strip()]
        if not input_dirs:
            st.error("Please provide at least one input directory.")
        else:
            combined_out = f"local_datasets/{output_name}/1_raw_combined"
            derep_out = f"local_datasets/{output_name}/2_dereplicated_genomes"
            status_container = st.empty()
            
            with st.spinner("Step 1: Physical Hashing and Cleaning..."):
                report_df = clean_and_rename_genomes(input_dirs, combined_out)
                
            if not report_df.empty:
                # Install/Check FastANI
                fastani_path = get_fastani_path()
                if not fastani_path:
                    with st.spinner("First time setup: Installing Pixi and FastANI..."):
                        try:
                            fastani_path = install_fastani(status_container)
                        except Exception as e:
                            st.error(f"Failed to install FastANI: {e}")
                            st.stop()
                            
                with st.spinner(f"Step 2: Running FastANI all-vs-all dereplication (Threshold: {derep_ani_thresh}% ANI, {derep_af_thresh}% AF)..."):
                    dup_report = run_fastani_dereplication(combined_out, derep_out, fastani_path, ani_threshold=derep_ani_thresh, af_threshold=derep_af_thresh, threads=threads, status_text=status_container)
                
                status_container.empty()
                st.success(f"Integration & Dereplication complete! Unique genomes saved to: `{os.path.abspath(derep_out)}`")
                
                col_s1, col_s2, col_s3, col_s4 = st.columns(4)
                total_files = report_df['Original_Count'].sum()
                physical_unique = len(report_df)
                physical_dups = total_files - physical_unique
                seq_dups = len(dup_report) if not dup_report.empty else 0
                final_unique = physical_unique - seq_dups
                
                col_s1.metric("Total Files Found", total_files)
                col_s2.metric("Physical Duplicates (MD5)", physical_dups)
                col_s3.metric("Sequence Duplicates (FastANI)", seq_dups)
                col_s4.metric("Final Unique Genomes", final_unique)
                
                if not dup_report.empty:
                    with st.expander("View FastANI Sequence Duplicates Removed"):
                        st.dataframe(dup_report)
                
            else:
                st.warning("No genome files (.fna, .fasta, .fa) found in the provided directories.")

    st.markdown("---")
    st.header("Step 3: Sequence Comparison vs Downloaded GTDB Data")
    st.markdown("Compare your dereplicated dataset directly against a folder of genomes you previously downloaded from GTDB (e.g., using the Custom Download tab).")
    
    comp_col1, comp_col2 = st.columns(2)
    with comp_col1:
        my_genomes_dir = st.text_input("Path to Your Unique Genomes Folder:", f"local_datasets/Bathyarchaeia_Combined/2_dereplicated_genomes")
    with comp_col2:
        gtdb_genomes_dir = st.text_input("Path to GTDB Downloaded Folder:", "ncbi_downloads/c__Bathyarchaeia_R232_all")
        
    comp_ani_thresh = st.slider("Comparison Match ANI Threshold (%)", min_value=90.0, max_value=100.0, value=95.0, step=0.5, help="Genomes with ANI >= this threshold will be considered the same species/matched.")
    comp_af_thresh = st.slider("Comparison Match AF Threshold (%)", min_value=10.0, max_value=100.0, value=60.0, step=1.0, help="Genomes must also have an Alignment Fraction (AF) >= this threshold to be considered matched.")
    
    if st.button("Run Sequence Comparison"):
        if not os.path.exists(my_genomes_dir) or not os.path.exists(gtdb_genomes_dir):
            st.error("One or both of the provided directory paths do not exist.")
        else:
            status_container = st.empty()
            fastani_path = get_fastani_path()
            if not fastani_path:
                with st.spinner("Installing FastANI..."):
                    fastani_path = install_fastani(status_container)
                    
            with st.spinner("Running FastANI Comparison... This may take a while depending on dataset size."):
                df_match = run_fastani_comparison(my_genomes_dir, gtdb_genomes_dir, fastani_path, ani_threshold=comp_ani_thresh, af_threshold=comp_af_thresh, threads=4, status_text=status_container)
                
            status_container.empty()
            
            if df_match is not None:
                my_total_files = set(os.listdir(my_genomes_dir))
                gtdb_total_files = set(os.listdir(gtdb_genomes_dir))
                
                # Exclude non-fasta files from count
                my_total = len([f for f in my_total_files if f.endswith(('.fna', '.fasta', '.fa'))])
                gtdb_total = len([f for f in gtdb_total_files if f.endswith(('.fna', '.fasta', '.fa'))])
                
                if df_match.empty:
                    matched_my = set()
                    matched_gtdb = set()
                else:
                    matched_my = set(df_match['query_name'])
                    matched_gtdb = set(df_match['ref_name'])
                    
                only_me = my_total - len(matched_my)
                only_gtdb = gtdb_total - len(matched_gtdb)
                
                st.subheader(f"Sequence Comparison Results (ANI ≥ {comp_ani_thresh}%, AF ≥ {comp_af_thresh}%)")
                
                st.markdown("##### 🧬 Your Local Dataset (My Genomes)")
                l1, l2, l3 = st.columns(3)
                l1.metric("Total Local Genomes", my_total)
                l2.metric("Matched Local Genomes", len(matched_my), help="Number of genomes in your dataset that found at least one match in GTDB.")
                l3.metric("Novel Local Genomes", only_me, help="Genomes in your dataset with NO match in GTDB.")
                
                st.markdown("##### 🌍 GTDB Downloaded Dataset")
                g1, g2, g3 = st.columns(3)
                g1.metric("Total GTDB Genomes", gtdb_total)
                g2.metric("Matched GTDB Genomes", len(matched_gtdb), help="Number of GTDB genomes that found at least one match in your dataset.")
                g3.metric("Missing GTDB Genomes", only_gtdb, help="Genomes in GTDB that have NO match in your dataset.")
                
                c_tab1, c_tab2, c_tab3 = st.tabs(["🤝 Matched Pairs (FastANI)", "🏠 Novel in My Dataset", "🌍 Missing (Only in GTDB)"])
                
                with c_tab1:
                    st.write(f"Found {len(df_match)} sequence match pairs crossing the threshold:")
                    st.dataframe(df_match[['query_name', 'ref_name', 'ani', 'af', 'matches', 'total']])
                    
                with c_tab2:
                    my_all = set([f for f in os.listdir(my_genomes_dir) if f.endswith(('.fna', '.fasta', '.fa'))])
                    novel_me = my_all - matched_my
                    st.write(f"These {len(novel_me)} genomes in your dataset have no match (ANI ≥ {comp_ani_thresh}%, AF ≥ {comp_af_thresh}%) in the GTDB folder:")
                    st.dataframe(pd.DataFrame({"Novel_Genome_File": list(novel_me)}))
                    
                with c_tab3:
                    gtdb_all = set([f for f in os.listdir(gtdb_genomes_dir) if f.endswith(('.fna', '.fasta', '.fa'))])
                    missing_gtdb = gtdb_all - matched_gtdb
                    st.write(f"These {len(missing_gtdb)} genomes in the GTDB folder have no match (ANI ≥ {comp_ani_thresh}%, AF ≥ {comp_af_thresh}%) in your dataset:")
                    st.dataframe(pd.DataFrame({"Missing_GTDB_File": list(missing_gtdb)}))
            else:
                st.error("Comparison failed or one of the directories contained no fasta files.")

