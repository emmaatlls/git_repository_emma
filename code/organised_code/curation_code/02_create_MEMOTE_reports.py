############################################################################################
# Imports
############################################################################################

import argparse
import cobra
import memote
import json
import pandas as pd
import re
import os

############################################################################################
# Paths
############################################################################################

DEFAULT_MODEL_DIR = os.path.join('/home/emma/Dokumente/thesis/Model_generation_curation/Curated_models/final')


def parse_args():
    parser = argparse.ArgumentParser(description='Generate MEMOTE reports for a directory of SBML models.')
    parser.add_argument('--model-dir', default=DEFAULT_MODEL_DIR,
                        help='Directory containing the SBML models to evaluate.')
    return parser.parse_args()

############################################################################################
# Configurations
############################################################################################

# Making sure the solver is set to cplex
cobra.Configuration().solver = "cplex"

############################################################################################
# Functions
############################################################################################

def create_reports(model_dir, report_dir):

    # Iterate through each model in the directory
    for model_file in os.listdir(model_dir):

        if model_file.endswith('.xml'):

            if f'{model_file[:-4]}_report.html' not in os.listdir(report_dir):

                # Read the model using cobrapy
                model = cobra.io.read_sbml_model(os.path.join(model_dir, model_file))

                # Test the model with MEMOTE
                model_results = memote.test_model(model, results=True)

                # Generate a html report
                model_report = memote.snapshot_report(model_results[1], html=True)

                # Save the report
                save_path = os.path.join(report_dir, f'{model_file[:-4]}_report.html')
                with open(save_path, 'w') as report_file:
                    report_file.write(model_report)

def extract_window_data(html_file):
    """Extract window.data dictionary from MEMOTE HTML report"""
    with open(html_file, 'r') as f:
        html_content = f.read()
	
    
    # Find window.data = {...}
    match = re.search(r'window\.data\s*=\s*({.*?});', html_content, re.DOTALL)
    if match:
        json_str = match.group(1)
        data_dict = json.loads(json_str)
        return data_dict
    else:
        raise ValueError("Could not find window.data in HTML file")

def summarize_reports(report_dir):
    data = []
    for file in os.listdir(report_dir):
        if file.endswith('.html'):
            html_file = os.path.join(report_dir, file)
            data_dict = extract_window_data(html_file)
            total_score = data_dict['score']['total_score']
            data.append({'Model': file[:-12], 'Total Score': total_score})
    summary = pd.DataFrame(data)
    summary.to_csv(os.path.join(report_dir, 'Summary_of_all_reports.csv'))

############################################################################################
# Main
############################################################################################
if __name__ == "__main__":
    args = parse_args()
    model_dir = os.path.abspath(args.model_dir)
    report_dir = os.path.join(model_dir, 'MEMOTE_reports')
    os.makedirs(report_dir, exist_ok=True)

    print(f"Checking target directory: {model_dir}")
    all_files = os.listdir(model_dir)
    print(f"Found total files in directory: {len(all_files)}")
    
    xml_files = [f for f in all_files if f.endswith(('.xml', '.sbml'))]
    print(f"Found {len(xml_files)} SBML/XML models to process.")

    if len(xml_files) > 0:
        print("Starting MEMOTE report generation...")
        create_reports(model_dir, report_dir)
        print("Starting report summarization...")
        summarize_reports(report_dir)
        print("Done! Check the MEMOTE_reports folder.")
    else:
        print("Execution stopped: No valid model files detected.")

