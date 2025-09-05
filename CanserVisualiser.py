import pandas as pd
import pymongo
from pymongo import MongoClient
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


MONGODB_CONNECTION_STRING = ""
MONGODB_DATABASE = "tcga_gene_expression"
MONGODB_COLLECTION = "patient_data"

print("Gene Expression Analysis & Visualization")
print("="*50)

# Initialize MongoDB client
try:
    mongo_client = MongoClient(MONGODB_CONNECTION_STRING)
    db = mongo_client[MONGODB_DATABASE]
    collection = db[MONGODB_COLLECTION]
    
    mongo_client.admin.command('ping')
    print(" MongoDB connection successful")
    
except Exception as e:
    print(f" MongoDB connection failed: {e}")
    exit(1)

def load_data_from_mongodb():
    """Load integrated data from MongoDB into pandas DataFrame"""
    print(" Loading integrated data from MongoDB...")
    
    # Get all documents with clinical data
    docs = list(collection.find({'clinical_data': {'$exists': True}}))
    
    if not docs:
        print(" No documents with clinical data found")
        return None
    
    print(f" Found {len(docs)} records with clinical data")
    
    # Convert to structured format
    rows = []
    
    for doc in docs:
        try:
            # Basic info
            row = {
                'patient_id': doc.get('patient_id'),
                'cancer_cohort': doc.get('cancer_cohort', '').replace('TCGA ', ''),
                'source_file': doc.get('source_file', '')
            }
            
            # Gene expression data
            gene_expr = doc.get('gene_expression', {})
            for gene, value in gene_expr.items():
                row[f'gene_{gene}'] = value
            
            # Clinical data
            clinical = doc.get('clinical_data', {})
            for key, value in clinical.items():
                row[f'clinical_{key}'] = value
            
            rows.append(row)
            
        except Exception as e:
            print(f"Error processing document: {e}")
            continue
    
    df = pd.DataFrame(rows)
    print(f" Loaded {len(df)} records into DataFrame")
    print(f" Columns: {len(df.columns)} total")
    
    return df

def analyze_survival_data(df):
    """Analyze survival outcomes"""
    print("\n Survival Analysis")
    print("-" * 30)
    
    # Overall Survival analysis
    if 'clinical_OS' in df.columns:
        os_counts = df['clinical_OS'].value_counts()
        print(f"Overall Survival (OS):")
        print(f"  Survived (1): {os_counts.get(1.0, 0)} patients")
        print(f"  Died (0): {os_counts.get(0.0, 0)} patients")
        print(f"  Missing: {df['clinical_OS'].isna().sum()} patients")
    
    # Disease Specific Survival analysis
    if 'clinical_DSS' in df.columns:
        dss_counts = df['clinical_DSS'].value_counts()
        print(f"\nDisease Specific Survival (DSS):")
        print(f"  Survived (1): {dss_counts.get(1.0, 0)} patients")
        print(f"  Died (0): {dss_counts.get(0.0, 0)} patients")
        print(f"  Missing: {df['clinical_DSS'].isna().sum()} patients")
    
    # Clinical stage analysis
    if 'clinical_clinical_stage' in df.columns:
        stage_counts = df['clinical_clinical_stage'].value_counts()
        print(f"\nClinical Stages:")
        for stage, count in stage_counts.head(10).items():
            print(f"  {stage}: {count} patients")

def analyze_gene_expression(df):
    """Analyze gene expression patterns"""
    print("\n Gene Expression Analysis")
    print("-" * 35)
    
    # Find gene expression columns
    gene_cols = [col for col in df.columns if col.startswith('gene_')]
    
    if not gene_cols:
        print(" No gene expression data found")
        return
    
    print(f" Analyzing {len(gene_cols)} genes")
    
    # Gene expression statistics
    for gene_col in gene_cols:
        gene_name = gene_col.replace('gene_', '')
        values = df[gene_col].dropna()
        
        if len(values) > 0:
            print(f"\n{gene_name}:")
            print(f"  Mean: {values.mean():.3f}")
            print(f"  Std: {values.std():.3f}")
            print(f"  Min: {values.min():.3f}")
            print(f"  Max: {values.max():.3f}")
            print(f"  Non-null samples: {len(values)}")

def create_visualizations(df):
    """Create comprehensive visualizations"""
    print("\n Creating Visualizations...")
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Gene expression columns
    gene_cols = [col for col in df.columns if col.startswith('gene_')]
    
    if len(gene_cols) == 0:
        print(" No gene expression data to visualize")
        return
    
    # 1. Gene Expression Distribution by Cancer Type
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('cGAS-STING Pathway Gene Expression Analysis', fontsize=16, fontweight='bold')
    
    # Select top 4 genes for visualization
    top_genes = gene_cols[:4]
    
    for i, gene_col in enumerate(top_genes):
        ax = axes[i//2, i%2]
        gene_name = gene_col.replace('gene_', '')
        
        # Box plot by cancer cohort
        cohort_data = []
        cohort_labels = []
        
        for cohort in df['cancer_cohort'].unique():
            if pd.notna(cohort):
                cohort_subset = df[df['cancer_cohort'] == cohort][gene_col].dropna()
                if len(cohort_subset) > 0:
                    cohort_data.append(cohort_subset)
                    cohort_labels.append(cohort[:10])  # Truncate long names
        
        if cohort_data:
            ax.boxplot(cohort_data, labels=cohort_labels)
            ax.set_title(f'{gene_name} Expression by Cancer Type')
            ax.set_ylabel('Expression Level')
            ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('gene_expression_by_cancer.png', dpi=300, bbox_inches='tight')
    print(" Saved: gene_expression_by_cancer.png")
    
    # 2. Survival Analysis Visualization
    if 'clinical_OS' in df.columns and len(gene_cols) > 0:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Gene Expression vs Overall Survival', fontsize=16, fontweight='bold')
        
        for i, gene_col in enumerate(top_genes):
            ax = axes[i//2, i%2]
            gene_name = gene_col.replace('gene_', '')
            
            # Create survival groups
            survival_data = df[[gene_col, 'clinical_OS']].dropna()
            
            if len(survival_data) > 0:
                survived = survival_data[survival_data['clinical_OS'] == 1.0][gene_col]
                died = survival_data[survival_data['clinical_OS'] == 0.0][gene_col]
                
                if len(survived) > 0 and len(died) > 0:
                    ax.hist([survived, died], bins=20, alpha=0.7, 
                           label=['Survived', 'Died'], color=['green', 'red'])
                    ax.set_title(f'{gene_name} Expression vs Survival')
                    ax.set_xlabel('Expression Level')
                    ax.set_ylabel('Number of Patients')
                    ax.legend()
                    
                    # Perform t-test
                    t_stat, p_value = stats.ttest_ind(survived, died)
                    ax.text(0.05, 0.95, f'p-value: {p_value:.4f}', 
                           transform=ax.transAxes, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig('gene_expression_vs_survival.png', dpi=300, bbox_inches='tight')
        print(" Saved: gene_expression_vs_survival.png")
    
    # 3. Correlation Heatmap
    if len(gene_cols) > 1:
        gene_data = df[gene_cols].select_dtypes(include=[np.number])
        
        if not gene_data.empty:
            # Calculate correlation matrix
            corr_matrix = gene_data.corr()
            
            plt.figure(figsize=(10, 8))
            mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
            sns.heatmap(corr_matrix, mask=mask, annot=True, cmap='coolwarm', center=0,
                       square=True, linewidths=0.5, cbar_kws={"shrink": .8})
            plt.title('Gene Expression Correlation Matrix', fontsize=14, fontweight='bold')
            
            # Clean up gene names for labels
            clean_labels = [col.replace('gene_', '') for col in corr_matrix.columns]
            plt.xticks(range(len(clean_labels)), clean_labels, rotation=45)
            plt.yticks(range(len(clean_labels)), clean_labels, rotation=0)
            
            plt.tight_layout()
            plt.savefig('gene_correlation_heatmap.png', dpi=300, bbox_inches='tight')
            print(" Saved: gene_correlation_heatmap.png")
    
    # 4. Clinical Stage Analysis
    if 'clinical_clinical_stage' in df.columns and len(gene_cols) > 0:
        # Filter out non-informative stages
        valid_stages = df[df['clinical_clinical_stage'].notna() & 
                         ~df['clinical_clinical_stage'].isin(['[Not Available]', '[Not Applicable]'])]
        
        if len(valid_stages) > 20:  # Only if we have enough data
            plt.figure(figsize=(12, 8))
            
            # Select first gene for stage analysis
            gene_col = gene_cols[0]
            gene_name = gene_col.replace('gene_', '')
            
            stage_order = ['Stage I', 'Stage II', 'Stage III', 'Stage IV']
            available_stages = [s for s in stage_order if s in valid_stages['clinical_clinical_stage'].values]
            
            if available_stages:
                sns.boxplot(data=valid_stages, x='clinical_clinical_stage', y=gene_col, 
                           order=available_stages)
                plt.title(f'{gene_name} Expression by Clinical Stage', fontsize=14, fontweight='bold')
                plt.xlabel('Clinical Stage')
                plt.ylabel('Expression Level')
                plt.xticks(rotation=45)
                
                plt.tight_layout()
                plt.savefig('gene_expression_by_stage.png', dpi=300, bbox_inches='tight')
                print("Saved: gene_expression_by_stage.png")
    
    plt.show()

def generate_summary_report(df):
    """Generate a comprehensive summary report"""
    print("\n COMPREHENSIVE ANALYSIS REPORT")
    print("=" * 50)
    
    # Dataset overview
    print(f"Dataset Overview:")
    print(f"  Total patients: {len(df)}")
    print(f"  Cancer cohorts: {df['cancer_cohort'].nunique()}")
    print(f"  Genes analyzed: {len([col for col in df.columns if col.startswith('gene_')])}")
    
    # Cohort breakdown
    print(f"\nCancer Cohort Distribution:")
    cohort_counts = df['cancer_cohort'].value_counts()
    for cohort, count in cohort_counts.items():
        print(f"  {cohort}: {count} patients")
    
    # Data completeness
    gene_cols = [col for col in df.columns if col.startswith('gene_')]
    if gene_cols:
        print(f"\nData Completeness:")
        for gene_col in gene_cols:
            gene_name = gene_col.replace('gene_', '')
            completeness = (1 - df[gene_col].isna().mean()) * 100
            print(f"  {gene_name}: {completeness:.1f}% complete")
    
    # Statistical insights
    if 'clinical_OS' in df.columns:
        print(f"\nSurvival Analysis Insights:")
        os_survival_rate = df['clinical_OS'].mean() * 100
        print(f"  Overall survival rate: {os_survival_rate:.1f}%")
    
    print(f"\n Analysis complete! Check the generated PNG files for visualizations.")

def main():
    """Main analysis function"""
    
    # Load data
    df = load_data_from_mongodb()
    if df is None:
        return
    
    # Perform analyses
    analyze_survival_data(df)
    analyze_gene_expression(df)
    
    # Create visualizations
    create_visualizations(df)
    
    # Generate summary report
    generate_summary_report(df)
    
    print(f"\n PROJECT COMPLETE!")
    print(f"Successfully analyzed gene expression data from {len(df)} patients")
    print(f" cGAS-STING pathway analysis finished")
    print(f" Check the generated visualizations for insights")

if __name__ == "__main__":
    main()