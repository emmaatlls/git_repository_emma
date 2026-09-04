#import argparse
import os
import subprocess

working_dir = '/home/emma/Dokumente/thesis/Model_generation_curation'
annotated_genomes_dir = os.path.join(working_dir, 'genome_files/faa_files')
gram_neg_dir = os.path.join(annotated_genomes_dir, 'gram_negative')
gram_pos_dir = os.path.join(annotated_genomes_dir, 'gram_positive')
gram_unknown_dir = os.path.join(annotated_genomes_dir, 'gram_unknown')
output_dir= os.path.join(working_dir, 'Draft_models/faa_models')
#DEFAULT_EXCLUDE_PATH = os.path.join(working_dir, 'exclude_IPPTandmodel.csv')
#exclude_path = os.path.join(working_dir, 'exclude_IPPTonly.csv')
#def parse_args():
#    parser = argparse.ArgumentParser(description='Carve draft models with an optional soft exclude file.')
 #   parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR,
 #                       help='Directory where carved models will be written.')
 #   parser.add_argument('--soft-exclude', default=DEFAULT_EXCLUDE_PATH,
 #                       help='Soft-exclude CSV file to pass to carve.')
 #   return parser.parse_args()


def carve_models(sequence_dir, universe, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for file in os.listdir(sequence_dir):
        if file.endswith('.faa'):
            sequence_path = os.path.join(sequence_dir, file)
            draft_path = os.path.join(output_dir, f'{file[:-4]}.xml')

            if os.path.exists(draft_path):
                print(f'Skipping {file}: Model already exists at {draft_path}')
                continue

            print(f'Starting carving process for {file}')
            subprocess.run([
                'conda', 'run', '-n', 'GEMS_creation',
                'carve', sequence_path,
                '--output', draft_path,
                '-u', universe
            ], check=True)
            print(f'Finished carving process for {file}')


def main():
    #args = parse_args()

    carve_models(gram_neg_dir, 'gramneg', output_dir)
    carve_models(gram_pos_dir, 'grampos', output_dir)

if __name__ == '__main__':
    main()
