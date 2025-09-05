from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from bs4 import BeautifulSoup
import time
import os
import requests
from urllib.parse import urljoin, urlparse
import boto3
from botocore.exceptions import ClientError
from io import BytesIO

AWS_ACCESS_KEY_ID = "NO"
AWS_SECRET_ACCESS_KEY = "NO"
S3_BUCKET_NAME = "gene-expression-data-mm"
AWS_REGION = "eu-north-1"

# Initialize S3 client
s3_client = boto3.client(
    's3',
    aws_access_key_id=AWS_ACCESS_KEY_ID,
    aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
    region_name=AWS_REGION
)

def download_and_upload_to_s3(url, bucket_name, object_name, also_save_local=False, local_dir=None):
    """
    Stream download directly to S3, optionally save local copy
    """
    try:
        # Handle relative URLs
        if not url.startswith('http'):
            url = urljoin('https://xenabrowser.net/', url)
        
        print(f"    Streaming from: {url}")
        
        # Stream the file
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()
        
        # Read into memory buffer
        file_buffer = BytesIO()
        file_size = 0
        
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                file_buffer.write(chunk)
                file_size += len(chunk)
        
        # Reset buffer position
        file_buffer.seek(0)
        
        # Upload to S3
        s3_client.upload_fileobj(file_buffer, bucket_name, object_name)
        print(f"    ✓ Streamed to S3: s3://{bucket_name}/{object_name} ({file_size:,} bytes)")
        
        # Optionally save local copy
        if also_save_local and local_dir:
            local_path = os.path.join(local_dir, os.path.basename(object_name))
            file_buffer.seek(0)  # Reset buffer
            
            with open(local_path, 'wb') as f:
                f.write(file_buffer.read())
            print(f"    ✓ Local copy: {local_path}")
        
        return True
        
    except Exception as e:
        print(f"    ✗ Error streaming {url}: {str(e)}")
        return False

def download_file_local_only(url, filepath):
    """Traditional download - keep for backup"""
    try:
        if not url.startswith('http'):
            url = urljoin('https://xenabrowser.net/', url)
            
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
            return True
        else:
            return False
            
    except Exception as e:
        print(f"    Error downloading {url}: {str(e)}")
        return False

# Setup Chrome options
options = Options()
options.headless = True
options.add_argument('--no-sandbox')
options.add_argument('--disable-dev-shm-usage')

service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=options)

# Create directory for optional local files
download_dir = 'gene_expression_data'
if not os.path.exists(download_dir):
    os.makedirs(download_dir)

urlbase = 'https://xenabrowser.net/datapages/'
urlOne = urlbase + '?hub=https://tcga.xenahubs.net:443'

print("Starting STREAMING gene expression scraper...")
print(f"S3 Bucket: {S3_BUCKET_NAME}")
print("Mode: Stream directly to S3 (no local storage required)")

driver.get(urlOne)
time.sleep(3)

html = driver.page_source
soup = BeautifulSoup(html, 'html.parser')

uploaded_files = []
processed_cohorts = []

ul_element = soup.find('ul')

if ul_element:
    li_elements = ul_element.find_all('li')
    print(f"Found {len(li_elements)} cancer cohorts to process...")
    
    for i, li in enumerate(li_elements):
        a_tag = li.find('a')
        if a_tag:
            cohort_name = a_tag.get_text(strip=True)
            link = a_tag['href']
            full_link = urlbase + link
            
            print(f"\n[{i+1}/{len(li_elements)}] Processing cohort: {cohort_name}")
            
            try:
                driver.get(full_link)
                time.sleep(3)
                
                new_page_html = driver.page_source
                new_page_soup = BeautifulSoup(new_page_html, 'html.parser')
                
                gene_expression_found = False
                all_links = new_page_soup.find_all('a', href=True)
                
                for link_tag in all_links:
                    link_text = link_tag.get_text(strip=True).lower()
                    
                    if ('illuminahiseq' in link_text and 
                        'pancan' in link_text and 
                        'normalized' in link_text):
                        
                        gene_expression_link = link_tag['href']
                        gene_expression_full_link = urlbase + gene_expression_link
                        
                        print(f"  Found gene expression data")
                        gene_expression_found = True
                        
                        driver.get(gene_expression_full_link)
                        time.sleep(3)
                        
                        gene_page_html = driver.page_source
                        gene_page_soup = BeautifulSoup(gene_page_html, 'html.parser')
                        
                        download_found = False
                        
                        # Look for download links
                        span_tags = gene_page_soup.find_all('span')
                        for span in span_tags:
                            if span.get_text() and 'download' in span.get_text().lower():
                                next_a = span.find_next('a', href=True)
                                if next_a:
                                    download_url = next_a['href']
                                    
                                    # Create S3 object name
                                    filename = f"{cohort_name.replace(' ', '_').replace('(', '').replace(')', '')}_gene_expression.tsv"
                                    s3_object_name = f"gene_expression/{filename}"
                                    
                                    print(f"  Streaming to S3...")
                                    
                                    # Stream directly to S3 (with optional local backup)
                                    if download_and_upload_to_s3(
                                        download_url, 
                                        S3_BUCKET_NAME, 
                                        s3_object_name,
                                        also_save_local=True,  # Change to False if you don't want local copies
                                        local_dir=download_dir
                                    ):
                                        uploaded_files.append({
                                            'cohort': cohort_name,
                                            'filename': filename,
                                            's3_path': f"s3://{S3_BUCKET_NAME}/{s3_object_name}",
                                            'url': download_url
                                        })
                                    
                                    download_found = True
                                    break
                        
                        # Fallback method
                        if not download_found:
                            download_links = gene_page_soup.find_all('a', href=True)
                            for dl_link in download_links:
                                href = dl_link['href']
                                if href.endswith('.tsv') or href.endswith('.txt') or 'download' in href.lower():
                                    download_url = href if href.startswith('http') else urljoin(gene_expression_full_link, href)
                                    
                                    filename = f"{cohort_name.replace(' ', '_').replace('(', '').replace(')', '')}_gene_expression.tsv"
                                    s3_object_name = f"gene_expression/{filename}"
                                    
                                    if download_and_upload_to_s3(
                                        download_url, 
                                        S3_BUCKET_NAME, 
                                        s3_object_name,
                                        also_save_local=True,
                                        local_dir=download_dir
                                    ):
                                        uploaded_files.append({
                                            'cohort': cohort_name,
                                            'filename': filename,
                                            's3_path': f"s3://{S3_BUCKET_NAME}/{s3_object_name}",
                                            'url': download_url
                                        })
                                        download_found = True
                                        break
                        
                        if not download_found:
                            print(f"  ⚠ No download link found")
                        
                        break
                
                if not gene_expression_found:
                    print(f"  ⚠ No gene expression data found")
                
                processed_cohorts.append({
                    'cohort': cohort_name,
                    'has_gene_expression': gene_expression_found
                })
                
                driver.back()
                time.sleep(2)
                
            except Exception as e:
                print(f"  ✗ Error processing {cohort_name}: {str(e)}")
                continue

# Summary
print("\n" + "="*60)
print("STREAMING UPLOAD SUMMARY")
print("="*60)
print(f"Total cohorts processed: {len(processed_cohorts)}")
print(f"Files successfully streamed to S3: {len(uploaded_files)}")

if uploaded_files:
    print(f"\nUploaded to S3 bucket '{S3_BUCKET_NAME}':")
    for file_info in uploaded_files:
        print(f"  - {file_info['cohort']}")
        print(f"    S3: {file_info['s3_path']}")

driver.quit()

print("\n Phase 2 Complete: All TSV files streamed to S3!")