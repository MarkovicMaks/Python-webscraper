from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from bs4 import BeautifulSoup
import time
import os
import requests
from urllib.parse import urljoin, urlparse

def download_file(url, filepath):
    """Download a file from URL and save it to filepath"""
    try:
        # Handle relative URLs
        if not url.startswith('http'):
            url = urljoin('https://xenabrowser.net/', url)
            
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        # Verify file was created and has content
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

# Auto-download and setup ChromeDriver - NO PATH NEEDED!
service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=options)

# Create directory for downloaded files
download_dir = 'gene_expression_data'
if not os.path.exists(download_dir):
    os.makedirs(download_dir)

urlbase = 'https://xenabrowser.net/datapages/'
urlOne = urlbase + '?hub=https://tcga.xenahubs.net:443'

print("Starting gene expression data scraper...")
print(f"Visiting main page: {urlOne}")

driver.get(urlOne)
time.sleep(3)

html = driver.page_source
soup = BeautifulSoup(html, 'html.parser')

# Find all cancer cohorts
ul_element = soup.find('ul')
downloaded_files = []
processed_cohorts = []

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
            print(f"Visiting: {full_link}")
            
            try:
                driver.get(full_link)
                time.sleep(3)
                
                new_page_html = driver.page_source
                new_page_soup = BeautifulSoup(new_page_html, 'html.parser')
                
                # Look for IlluminaHiSeq pancan normalized in gene expression section
                gene_expression_found = False
                
                # Find all links that contain gene expression data
                all_links = new_page_soup.find_all('a', href=True)
                
                for link_tag in all_links:
                    link_text = link_tag.get_text(strip=True).lower()
                    
                    # Look for IlluminaHiSeq pancan normalized files
                    if ('illuminahiseq' in link_text and 
                        'pancan' in link_text and 
                        'normalized' in link_text):
                        
                        gene_expression_link = link_tag['href']
                        gene_expression_full_link = urlbase + gene_expression_link
                        
                        print(f"  Found gene expression data: {gene_expression_full_link}")
                        gene_expression_found = True
                        
                        # Navigate to the gene expression page
                        driver.get(gene_expression_full_link)
                        time.sleep(3)
                        
                        gene_page_html = driver.page_source
                        gene_page_soup = BeautifulSoup(gene_page_html, 'html.parser')
                        
                        # Look for download link
                        download_found = False
                        
                        # Method 1: Look for span with "download" text
                        span_tags = gene_page_soup.find_all('span')
                        for span in span_tags:
                            if span.get_text() and 'download' in span.get_text().lower():
                                # Find next anchor tag
                                next_a = span.find_next('a', href=True)
                                if next_a:
                                    download_url = next_a['href']
                                    
                                    # Create filename based on cohort
                                    filename = f"{cohort_name.replace(' ', '_').replace('(', '').replace(')', '')}_gene_expression.tsv"
                                    filepath = os.path.join(download_dir, filename)
                                    
                                    print(f"  Downloading: {download_url}")
                                    print(f"  Saving as: {filename}")
                                    
                                    # Download the file
                                    if download_file(download_url, filepath):
                                        downloaded_files.append({
                                            'cohort': cohort_name,
                                            'filename': filename,
                                            'filepath': filepath,
                                            'url': download_url
                                        })
                                        print(f"  ✓ Successfully downloaded")
                                    else:
                                        print(f"  ✗ Failed to download")
                                    
                                    download_found = True
                                    break
                        
                        # Method 2: Look for direct download links if Method 1 fails
                        if not download_found:
                            download_links = gene_page_soup.find_all('a', href=True)
                            for dl_link in download_links:
                                href = dl_link['href']
                                if href.endswith('.tsv') or href.endswith('.txt') or 'download' in href.lower():
                                    download_url = href if href.startswith('http') else urljoin(gene_expression_full_link, href)
                                    
                                    filename = f"{cohort_name.replace(' ', '_').replace('(', '').replace(')', '')}_gene_expression.tsv"
                                    filepath = os.path.join(download_dir, filename)
                                    
                                    print(f"  Downloading (method 2): {download_url}")
                                    print(f"  Saving as: {filename}")
                                    
                                    if download_file(download_url, filepath):
                                        downloaded_files.append({
                                            'cohort': cohort_name,
                                            'filename': filename,
                                            'filepath': filepath,
                                            'url': download_url
                                        })
                                        print(f"  ✓ Successfully downloaded")
                                        download_found = True
                                        break
                                    else:
                                        print(f"  ✗ Failed to download")
                        
                        if not download_found:
                            print(f"  ⚠ No download link found for gene expression data")
                        
                        break  # Exit after finding first valid gene expression dataset
                
                if not gene_expression_found:
                    print(f"  ⚠ No IlluminaHiSeq pancan normalized gene expression data found for {cohort_name}")
                
                processed_cohorts.append({
                    'cohort': cohort_name,
                    'has_gene_expression': gene_expression_found
                })
                
                # Go back to cohort list
                driver.back()
                time.sleep(2)
                
            except Exception as e:
                print(f"  ✗ Error processing {cohort_name}: {str(e)}")
                continue

else:
    print("The <ul> element was not found.")

# Summary
print("\n" + "="*60)
print("SCRAPING SUMMARY")
print("="*60)
print(f"Total cohorts processed: {len(processed_cohorts)}")
print(f"Cohorts with gene expression data: {len([c for c in processed_cohorts if c['has_gene_expression']])}")
print(f"Files successfully downloaded: {len(downloaded_files)}")

if downloaded_files:
    print("\nDownloaded files:")
    for file_info in downloaded_files:
        print(f"  - {file_info['cohort']}: {file_info['filename']}")
        
print(f"\nAll files saved in: {os.path.abspath(download_dir)}")

driver.quit()