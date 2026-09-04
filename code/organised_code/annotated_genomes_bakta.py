############################################################################################
# Imports
############################################################################################

import os 
import subprocess

############################################################################################
# Functions
############################################################################################

def check_create_paths(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

############################################################################################
# Paths
############################################################################################

working_dir = "/home/emma/Dokumente/thesis/Model_generation_curation/genome_files"
annotated_neg_genomes_dir = os.path.join(working_dir, 'gram_negative')
annotated_pos_genomes_dir = os.path.join(working_dir, 'gram_positive')
annotated_protein_files_dir = os.path.join(working_dir, "faa_files")

# Ensure output directory exists
check_create_paths(annotated_protein_files_dir)
db_dir = os.path.join(working_dir, "bakta_db/")
target_directories = [annotated_neg_genomes_dir, annotated_pos_genomes_dir]

############################################################################################
# Processing
############################################################################################

for target_dir in target_directories:
    if not os.path.exists(target_dir):
        continue
        
    for file_name in os.listdir(target_dir):
        if file_name.endswith('.fa'):
            sample_contig = os.path.join(target_dir, file_name)
            
            sample_name = os.path.splitext(file_name)[0]
            
            sample_annotation = os.path.join(annotated_protein_files_dir, f'{sample_name}_bakta')
            
            print(f"Running Bakta annotation for: {sample_name}...")
            
            try: #use bakta to turn .fa into .faa files
                subprocess.run([
                    'conda', 'run', '-n', 'gapfill',
                    'bakta', 
                    '--threads', '8', 
                    '--prefix', sample_name,
                    "--skip-sorf", 
                    '--output', sample_annotation, 
                    '--db', db_dir,
                    sample_contig
                ], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error annotating {sample_name}: {e}")