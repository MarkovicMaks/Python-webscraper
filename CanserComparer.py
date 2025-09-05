import pandas as pd
import pymongo
from pymongo import MongoClient
import numpy as np


MONGODB_CONNECTION_STRING = "" 
MONGODB_DATABASE = "tcga_gene_expression"
MONGODB_COLLECTION = "patient_data"

CLINICAL_DATA_PATH = "TCGA_clinical_survival_data.tsv" 

print("Clinical Data Integration")
print("="*40)

try:
    mongo_client = MongoClient(MONGODB_CONNECTION_STRING)
    db = mongo_client[MONGODB_DATABASE]
    collection = db[MONGODB_COLLECTION]
    
    # Test connection
    mongo_client.admin.command('ping')
    print("MongoDB connection successful")
    
except Exception as e:
    print(f"MongoDB connection failed: {e}")
    exit(1)

def load_clinical_data():
    """Load and process clinical survival data"""
    try:
        print(f"Loading clinical data from: {CLINICAL_DATA_PATH}")
        
        # Read the clinical data
        clinical_df = pd.read_csv(CLINICAL_DATA_PATH, sep='\t', low_memory=False)
        
        print(f" Clinical data shape: {clinical_df.shape}")
        print(f" Columns: {list(clinical_df.columns[:10])}...")
        
        # Select only the columns we need
        required_columns = [
            'bcr_patient_barcode',  # Patient ID for joining
            'DSS',                  # Disease Specific Survival
            'OS',                   # Overall Survival  
            'clinical_stage'        # Clinical stage
        ]
        
        # Check which columns exist
        existing_columns = []
        for col in required_columns:
            if col in clinical_df.columns:
                existing_columns.append(col)
            else:
                print(f"  Column '{col}' not found in clinical data")
        
        if not existing_columns:
            print(" No required columns found in clinical data")
            return None
            
        # Select only existing columns
        clinical_subset = clinical_df[existing_columns].copy()
        
        # Clean up patient IDs - remove any extra characters
        if 'bcr_patient_barcode' in clinical_subset.columns:
            clinical_subset['bcr_patient_barcode'] = clinical_subset['bcr_patient_barcode'].astype(str).str.strip()
        
        print(f" Loaded {len(clinical_subset)} clinical records")
        print(f" Available columns: {existing_columns}")
        
        return clinical_subset
        
    except Exception as e:
        print(f"Error loading clinical data: {e}")
        return None

def extract_patient_id_from_tcga(tcga_id):
    """Extract patient barcode from TCGA sample ID"""
    try:
        # TCGA IDs are like: TCGA-AB-2872-03
        # Patient barcode is: TCGA-AB-2872
        parts = str(tcga_id).split('-')
        if len(parts) >= 3:
            return f"{parts[0]}-{parts[1]}-{parts[2]}"
        return str(tcga_id)
    except:
        return str(tcga_id)

def integrate_clinical_data(clinical_df):
    """Integrate clinical data with existing gene expression data in MongoDB"""
    
    print("\n Starting clinical data integration...")
    
    # Get all documents from MongoDB
    print("Fetching gene expression data from MongoDB...")
    gene_expression_docs = list(collection.find())
    
    print(f"  Found {len(gene_expression_docs)} gene expression records")
    
    if not gene_expression_docs:
        print(" No gene expression data found in MongoDB")
        return
    
    # Convert clinical data to dictionary for faster lookup
    clinical_dict = {}
    if 'bcr_patient_barcode' in clinical_df.columns:
        for _, row in clinical_df.iterrows():
            patient_id = str(row['bcr_patient_barcode']).strip()
            clinical_dict[patient_id] = row.to_dict()
    
    print(f"  Clinical data available for {len(clinical_dict)} patients")
    
    # Process each gene expression document
    updated_count = 0
    matched_count = 0
    
    for doc in gene_expression_docs:
        try:
            # Extract patient barcode from TCGA ID
            tcga_patient_id = doc.get('patient_id', '')
            patient_barcode = extract_patient_id_from_tcga(tcga_patient_id)
            
            # Look for clinical data
            clinical_data = clinical_dict.get(patient_barcode)
            
            if clinical_data:
                matched_count += 1
                
                # Prepare clinical data for update (exclude patient ID)
                clinical_update = {}
                for key, value in clinical_data.items():
                    if key != 'bcr_patient_barcode':
                        # Handle NaN values
                        if pd.isna(value):
                            clinical_update[key] = None
                        else:
                            clinical_update[key] = value
                
                # Update MongoDB document
                collection.update_one(
                    {'_id': doc['_id']},
                    {'$set': {'clinical_data': clinical_update}}
                )
                
                updated_count += 1
                
                if updated_count % 100 == 0:
                    print(f"   💾 Updated {updated_count} records...")
            
        except Exception as e:
            print(f"   ⚠️ Error processing document {doc.get('_id')}: {e}")
            continue
    
    print(f"\n Clinical data integration complete!")
    print(f"   Total gene expression records: {len(gene_expression_docs)}")
    print(f"   Records with matching clinical data: {matched_count}")
    print(f"   Records updated: {updated_count}")
    
    return updated_count

def show_sample_integrated_data():
    """Display sample of integrated data"""
    print("\n📋 Sample integrated data:")
    
    # Find documents with clinical data
    sample_docs = list(collection.find(
        {'clinical_data': {'$exists': True}},
        limit=3
    ))
    
    for i, doc in enumerate(sample_docs, 1):
        print(f"\nSample {i}:")
        print(f"  Patient ID: {doc.get('patient_id')}")
        print(f"  Cancer Cohort: {doc.get('cancer_cohort')}")
        
        # Show gene expression data
        gene_expr = doc.get('gene_expression', {})
        print(f"  Gene Expression: {len(gene_expr)} genes")
        for gene, value in list(gene_expr.items())[:3]:
            print(f"    {gene}: {value}")
        
        # Show clinical data
        clinical = doc.get('clinical_data', {})
        print(f"  Clinical Data:")
        for key, value in clinical.items():
            print(f"    {key}: {value}")

def main():
    """Main integration function"""
    
    # Load clinical data
    clinical_df = load_clinical_data()
    if clinical_df is None:
        return
    
    # Integrate with gene expression data
    updated_count = integrate_clinical_data(clinical_df)
    
    if updated_count > 0:
        # Show sample results
        show_sample_integrated_data()
        
        print(f"\n Integration successful!")
        print(f"{updated_count} records now have both gene expression and clinical data")
        print(f"Ready for visualization and analysis!")
    else:
        print(" No records were updated")

if __name__ == "__main__":
    main()