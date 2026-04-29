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

def clean_and_rename_genomes(input_dirs, output_main_dir, prefix="GENOME"):
    """
    Clean, rename and organize genomes from multiple directories.
    Uses MD5 hashing to skip physical duplicates.
    Prioritizes files that look like NCBI accessions (GCA/GCF).
    """
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

    genome_counter = 1
    for f_hash, file_list in hash_to_files.items():
        # Heuristic: Pick the "best" filename among duplicates
        # Priority: 1. Contains GCA/GCF, 2. Shortest name (often cleaner)
        best_file = file_list[0]
        for f in file_list:
            if re.search(r'GC[AF]_\d+\.\d+', f['name']):
                best_file = f
                break
        
        new_filename = f"{prefix}_{genome_counter:04d}.fasta"
        new_path = os.path.join(output_main_dir, new_filename)
        
        shutil.copy2(best_file['path'], new_path)
        
        report.append({
            "Selected_Source": best_file['name'],
            "Original_Count": len(file_list),
            "New_Name": new_filename,
            "Other_Names": ", ".join([f['name'] for f in file_list if f['name'] != best_file['name']])
        })
        genome_counter += 1
                
    return pd.DataFrame(report)

def extract_accession_from_filename(name):
    """Extract GCA/GCF accession from string."""
    match = re.search(r'(GC[AF]_\d+\.\d+)', name)
    return match.group(1) if match else None

def run_drep(input_dir, output_dir, comp_threshold=50, cont_threshold=10):
    """Placeholder for dRep execution."""
    # This would normally be: subprocess.run(["dRep", "dereplicate", output_dir, "-g", input_dir + "/*.fasta", ...])
    # Since dRep is not in the environment, we'll provide the command for the user
    cmd = f"dRep dereplicate {output_dir} -g {input_dir}/*.fasta --completeness {comp_threshold} --contamination {cont_threshold}"
    return cmd

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
        st.write(f"Found **{len(df_res)}** genomes matching `{search_single}`.")
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
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Genomes in R" + str(v1), len(df_v1))
            col2.metric("Genomes in R" + str(v2), len(df_v2))
            col3.metric("Net Change", len(df_v2) - len(df_v1))
            
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
    st.header("Local Dataset Updater & Dereplication")
    st.markdown("""
    This module helps you integrate genomes from multiple local folders, perform basic cleaning/renaming, 
    and provides commands for dereplication (dRep).
    """)
    
    with st.expander("Step 1: Configure Input Directories", expanded=True):
        input_dirs_raw = st.text_area("Enter paths to local genome folders (one per line):", 
                                     help="Example: /Users/user/genomes/bathy_v1\n/Users/user/genomes/bathy_v2")
        output_name = st.text_input("New Dataset Name:", "Bathyarchaeia_Combined")
        genome_prefix = st.text_input("Genome ID Prefix:", "BATHY")
        
    if st.button("Analyze & Integrate Datasets"):
        input_dirs = [d.strip() for d in input_dirs_raw.split('\n') if d.strip()]
        if not input_dirs:
            st.error("Please provide at least one input directory.")
        else:
            combined_out = f"local_datasets/{output_name}/cleaned_genomes"
            
            with st.spinner("Processing genomes (Hashing, Cleaning, Renaming)..."):
                report_df = clean_and_rename_genomes(input_dirs, combined_out, prefix=genome_prefix)
                
            if not report_df.empty:
                st.success(f"Integration complete! Genomes saved to: `{os.path.abspath(combined_out)}`")
                
                # Display Summary Statistics
                col_s1, col_s2, col_s3 = st.columns(3)
                total_found = len(report_df)
                success_count = len(report_df[report_df['Status'] == 'Success'])
                dup_count = len(report_df[report_df['Status'].str.contains('Duplicate')])
                
                col_s1.metric("Total Files Found", total_found)
                col_s2.metric("Unique Genomes", success_count)
                col_s3.metric("Exact Duplicates", dup_count)
                
                st.subheader("Integration Report")
                st.dataframe(report_df)
                
                # Step 2: dRep Integration
                st.header("Step 2: Dereplication (dRep)")
                st.info("""
                dRep is used to identify representative genomes based on Average Nucleotide Identity (ANI).
                Since dRep requires a complex environment (MASH, MUMmer, etc.), please run the following command in your terminal:
                """)
                
                drep_out = f"local_datasets/{output_name}/drep_results"
                drep_cmd = run_drep(os.path.abspath(combined_out), os.path.abspath(drep_out))
                
                st.code(drep_cmd, language="bash")
                st.markdown(f"**Output directory will be:** `{os.path.abspath(drep_out)}`")
                
                # --- New Dataset vs GTDB Comparison Section ---
                st.header("Step 3: Compare with GTDB Taxon")
                st.markdown("Compare your integrated dataset with a specific GTDB release and taxon.")
                
                comp_col1, comp_col2 = st.columns(2)
                with comp_col1:
                    gtdb_v = st.selectbox("Select GTDB Version to compare", versions, index=len(versions)-1, key="upd_gtdb_v")
                with comp_col2:
                    gtdb_t = st.text_input("Enter GTDB Taxon (e.g., c__Bathyarchaeia)", "c__Bathyarchaeia", key="upd_gtdb_t")
                
                if st.button("Compare Dataset vs GTDB"):
                    # 1. Get GTDB Accessions for that taxon
                    df_gtdb = df[(df['Version'] == gtdb_v) & (df['Taxonomy'].str.contains(gtdb_t, na=False))]
                    gtdb_accs = set([g[3:] if g.startswith(("RS_", "GB_")) else g for g in df_gtdb['Genome_ID']])
                    
                    # 2. Get My Dataset Accessions (from report_df)
                    # We look at all filenames in the original folders that were mapped to this dataset
                    my_accs = set()
                    for idx, row in report_df.iterrows():
                        # Check selected source and other names for accessions
                        all_names = [row['Selected_Source']] + (row['Other_Names'].split(", ") if row['Other_Names'] else [])
                        for name in all_names:
                            acc = extract_accession_from_filename(name)
                            if acc:
                                my_accs.add(acc)
                    
                    if not my_accs:
                        st.warning("No NCBI Accessions (GCA/GCF) detected in your local filenames. Comparison is only possible if filenames contain accessions.")
                    else:
                        common = my_accs.intersection(gtdb_accs)
                        only_me = my_accs - gtdb_accs
                        only_gtdb = gtdb_accs - my_accs
                        
                        st.subheader(f"Composition Comparison: My Dataset vs GTDB R{gtdb_v} ({gtdb_t})")
                        
                        m1, m2, m3 = st.columns(3)
                        m1.metric("Common (Overlap)", len(common))
                        m2.metric("Only in My Dataset", len(only_me))
                        m3.metric("Only in GTDB", len(only_gtdb))
                        
                        c_tab1, c_tab2, c_tab3 = st.tabs(["🤝 Common", "🏠 Only Me", "🌍 Only GTDB"])
                        
                        with c_tab1:
                            st.write(f"These {len(common)} genomes are present in both sets.")
                            st.dataframe(pd.DataFrame({"Accession": list(common)}))
                            
                        with c_tab2:
                            st.write(f"These {len(only_me)} genomes are in your local folders but NOT in the GTDB {gtdb_t} group.")
                            st.dataframe(pd.DataFrame({"Accession": list(only_me)}))
                            
                        with c_tab3:
                            st.write(f"These {len(only_gtdb)} genomes are in GTDB {gtdb_t} but NOT in your local folders.")
                            st.dataframe(pd.DataFrame({"Accession": list(only_gtdb)}))
            else:
                st.warning("No genome files (.fna, .fasta, .fa) found in the provided directories.")

