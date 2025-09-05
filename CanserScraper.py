from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import time

#ne znam zašto ipak prikazuje ui
options = Options()
options.headless = True

#tu promjeni path do svojeg chrome drivera
chrome_driver_path = 'C:/Users/Korisnik/Desktop/chromedriver-win64/chromedriver.exe'

service = Service(chrome_driver_path)
driver = webdriver.Chrome(service=service, options=options)

urlbase = 'https://xenabrowser.net/datapages/'
urlOne = urlbase + '?hub=https://tcga.xenahubs.net:443'
driver.get(urlOne)
time.sleep(3) 

html = driver.page_source
soup = BeautifulSoup(html, 'html.parser')

ul_element = soup.find('ul')
if ul_element:
    li_elements = ul_element.find_all('li')

    for li in li_elements:
        a_tag = li.find('a')
        if a_tag:
            link = a_tag['href']
            full_link = urlbase + link
            print(f"Visiting: {full_link}")

            driver.get(full_link)
            time.sleep(3) 
            
            new_page_html = driver.page_source
            new_page_soup = BeautifulSoup(new_page_html, 'html.parser')

            
            pancancake_link_tag = new_page_soup.find('a', string=lambda text: text and text.endswith('pancan normalized'))
            if pancancake_link_tag:
                pancancake_link = pancancake_link_tag['href']
                pancancake_link_Full = urlbase + pancancake_link
                print(f"Found pancancake normalized link: {pancancake_link_Full}")

                driver.get(pancancake_link_Full)
                time.sleep(3)
                
                pancancake_page_html = driver.page_source
                pancancake_page_soup = BeautifulSoup(pancancake_page_html, 'html.parser')

                span_tag = pancancake_page_soup.find('span', string=lambda text: text and 'download' in text.lower())

                if span_tag:
                    # zaš sam sve nije moglo biti u jednom elementu
                    a_tag = span_tag.find_next('a', href=True)
                    if a_tag:
                        download_link = a_tag['href']
                        driver.get(download_link)
                        time.sleep(2)
                        driver.back()
                    else:
                        print("No <a> tag found near the 'download' span.")
                else:
                    print("No 'download' span found.")
                
            else:
                print("No pancancake normalized link found on this page.")
            driver.back()
            time.sleep(2)
else:
    print("The <ul> element was not found.")


driver.quit()
