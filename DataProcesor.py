import boto3
import pandas as pd
import pymongo
from pymongo import MongoClient
from io import StringIO
import numpy as np
import time
import gzip


AWS_ACCESS_KEY_ID = ""
AWS_SECRET_ACCESS_KEY = ""
S3_BUCKET_NAME = "gene-expression-data-mm"
AWS_REGION = "eu-north-1"

MONGODB_CONNECTION_STRING = "" 
MONGODB_DATABASE = "tcga_gene_expression"
MONGODB_COLLECTION = "patient_data"


CGAS_STING_GENES = [
    'C6orf150',  
    'CXCL10', 
    'TMEM173',   
    'CXCL9',
    'CXCL11',
    'NFKB1',
    'IKBKE',
    'IRF3',
    'TREX1',
    'ATM',
    'IL6',
    'IL8',      
    'CXCL8'    
]

print("🧬 TCGA Gene Expression Data Processor")
print("="*50)

s3_client = boto3.client(
    's3',
    aws_access_key_id=AWS_ACCESS_KEY_ID,
    aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
    region_name=AWS_REGION
)

try:
    mongo_client = MongoClient(MONGODB_CONNECTION_STRING)
    db = mongo_client[MONGODB_DATABASE]
    collection = db[MONGODB_COLLECTION]
    
    # Test connection
    mongo_client.admin.command('ping')
    print(" MongoDB connection successful")
    
except Exception as e:
    print(f" MongoDB connection failed: {e}")
    exit(1)

def list_s3_files():
    """List all TSV files in the S3 bucket"""
    try:
        response = s3_client.list_objects_v2(
            Bucket=S3_BUCKET_NAME,
            Prefix='gene_expression/'
        )
        
        files = []
        if 'Contents' in response:
            for obj in response['Contents']:
                if obj['Key'].endswith('.tsv'):
                    files.append(obj['Key'])
        
        return files
    
    except Exception as e:
        print(f" Error listing S3 files: {e}")
        return []

def read_tsv_from_s3(s3_key):
    """Read TSV file from S3 into pandas DataFrame (handles gzipped files)"""
    try:
        print(f" Reading: {s3_key}")
        
        # Download file content
        response = s3_client.get_object(Bucket=S3_BUCKET_NAME, Key=s3_key)
        content_bytes = response['Body'].read()
        
        # Check if file is gzipped (starts with magic bytes 0x1f 0x8b)
        if content_bytes[:2] == b'\x1f\x8b':
            print(f"   🗜️ File is gzipped, decompressing...")
            content = gzip.decompress(content_bytes).decode('utf-8')
        else:
            print(f"    File is plain text")
            content = content_bytes.decode('utf-8')
        
        # Read into pandas DataFrame
        df = pd.read_csv(StringIO(content), sep='\t', low_memory=False)
        
        print(f"   Shape: {df.shape}")
        print(f"   Columns: {list(df.columns[:5])}...")  # Show first 5 columns
        
        return df
        
    except Exception as e:
        print(f"Error reading {s3_key}: {e}")
        return None

def extract_cohort_from_filename(filename):
    """Extract cancer cohort name from filename"""
    # Remove path and extension, clean up
    base_name = filename.split('/')[-1].replace('_gene_expression.tsv', '')
    return base_name.replace('_', ' ')

def process_gene_expression_file(s3_key):
    """Process a single gene expression file"""
    
    # Read the TSV file
    df = read_tsv_from_s3(s3_key)
    if df is None:
        return 0
   
    cohort = extract_cohort_from_filename(s3_key)

    gene_name_col = df.columns[0]  # Usually 'sample'
    df_transposed = df.set_index(gene_name_col)
    
    print(f"   Data structure: {df_transposed.shape[0]} genes × {df_transposed.shape[1]} patients")
    
    available_genes = []
    gene_rows = {}
    
    for target_gene in CGAS_STING_GENES:
        matching_genes = [gene for gene in df_transposed.index if target_gene.upper() in str(gene).upper()]
        
        if matching_genes:
            gene_rows[target_gene] = matching_genes[0]
            available_genes.append(target_gene)
    
    print(f"   🎯 Found {len(available_genes)} target genes: {available_genes}")
    
    if not available_genes:
        print(f"   ⚠️ No target genes found in {s3_key}")
        print(f"   🔍 Sample gene names: {list(df_transposed.index[:10])}")
        return 0
    
    records_inserted = 0
    batch_records = []
    
    for patient_id in df_transposed.columns:
        try:
            if pd.isna(patient_id) or patient_id == '':
                continue
            
            gene_expression = {}
            for gene, gene_name in gene_rows.items():
                value = df_transposed.loc[gene_name, patient_id]
                
                if pd.isna(value):
                    gene_expression[gene] = None
                else:
                    try:
                        gene_expression[gene] = float(value)
                    except (ValueError, TypeError):
                        gene_expression[gene] = None
            
            # Create document for MongoDB
            document = {
                'patient_id': str(patient_id),
                'cancer_cohort': cohort,
                'gene_expression': gene_expression,
                'processing_date': time.time(),
                'source_file': s3_key
            }
            
            batch_records.append(document)
            
            if len(batch_records) >= 500:
                try:
                    collection.insert_many(batch_records)
                    records_inserted += len(batch_records)
                    print(f"   Inserted batch: {records_inserted} records")
                    batch_records = []
                except Exception as e:
                    print(f"   Batch insert error: {e}")
                    batch_records = []
        
        except Exception as e:
            print(f"   Error processing patient {patient_id}: {e}")
            continue
    
    if batch_records:
        try:
            collection.insert_many(batch_records)
            records_inserted += len(batch_records)
        except Exception as e:
            print(f"   Final batch insert error: {e}")

    print(f"  Processed {records_inserted} patients from {cohort}")
    return records_inserted

def create_mongodb_indexes():
    """Create indexes for better query performance"""
    try:
        collection.create_index("patient_id")
        collection.create_index("cancer_cohort")
        collection.create_index([("patient_id", 1), ("cancer_cohort", 1)])
        print(" MongoDB indexes created")
    except Exception as e:
        print(f" Error creating indexes: {e}")

def main():
    """Main processing function"""
    
    print("\n Step 1: Listing files in S3...")
    s3_files = list_s3_files()
    
    if not s3_files:
        print("❌ No TSV files found in S3")
        return
    
    print(f" Found {len(s3_files)} TSV files")
    for file in s3_files:
        print(f"   - {file}")
    
    print("\n Step 2: Setting up MongoDB...")
    create_mongodb_indexes()
    
    # Clear existing data (optional - remove if you want to keep existing data)
    print("Clearing existing data...")
    collection.delete_many({})
    
    print("\n Step 3: Processing files...")
    total_records = 0
    
    for i, s3_key in enumerate(s3_files, 1):
        print(f"\n[{i}/{len(s3_files)}] Processing: {s3_key}")
        
        try:
            records = process_gene_expression_file(s3_key)
            total_records += records
            
        except Exception as e:
            print(f" Failed to process {s3_key}: {e}")
            continue
    
    print("\n" + "="*50)
    print(" PROCESSING COMPLETE!")
    print("="*50)
    print(f" Total records inserted: {total_records:,}")
    print(f" Files processed: {len(s3_files)}")
    
    # Show sample data
    print("\n Sample documents in MongoDB:")
    sample_docs = list(collection.find().limit(3))
    for i, doc in enumerate(sample_docs, 1):
        print(f"\nSample {i}:")
        print(f"  Patient ID: {doc.get('patient_id')}")
        print(f"  Cohort: {doc.get('cancer_cohort')}")
        print(f"  Genes: {list(doc.get('gene_expression', {}).keys())}")
    
    print(f"\n Data ready for clinical data integration!")

if __name__ == "__main__":
    main()