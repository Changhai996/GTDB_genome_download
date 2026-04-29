import os
import re
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin

BASE_URL = "https://data.gtdb.ecogenomic.org/releases/"
DOWNLOAD_DIR = "gtdb_data"

if not os.path.exists(DOWNLOAD_DIR):
    os.makedirs(DOWNLOAD_DIR)

def get_links(url):
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        links = []
        for a in soup.find_all('a'):
            href = a.get('href')
            if href and not href.startswith('?') and href != '/' and not href.startswith('http'):
                links.append(href)
        return links
    except Exception as e:
        print(f"Error fetching {url}: {e}")
        return []

def crawl_releases():
    releases = get_links(BASE_URL)
    for rel in releases:
        if not rel.startswith('release'):
            continue
        rel_url = urljoin(BASE_URL, rel)
        versions = get_links(rel_url)
        for ver in versions:
            # Look for version folders like 202.0/
            if re.match(r'^\d+\.\d+/?$', ver):
                ver_url = urljoin(rel_url, ver)
                files = get_links(ver_url)
                for file in files:
                    # We are looking for ar*_taxonomy*.tsv or .tsv.gz
                    if re.match(r'^ar\d+_taxonomy_r\d+\.tsv(\.gz)?$', file):
                        file_url = urljoin(ver_url, file)
                        download_file(file_url, file)

def download_file(url, filename):
    filepath = os.path.join(DOWNLOAD_DIR, filename)
    if os.path.exists(filepath):
        print(f"Already downloaded: {filename}")
        return
    
    print(f"Downloading {filename} from {url} ...")
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(filepath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded {filename}")
    except Exception as e:
        print(f"Failed to download {filename}: {e}")

if __name__ == "__main__":
    crawl_releases()
