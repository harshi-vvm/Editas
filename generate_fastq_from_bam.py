import subprocess
import os
import pandas as pd
import jira_operations



def load_json_file(file_path):

    with open(file_path, 'r') as file:
        config = json.load(file)
    return config


def download_analysis_sheet(data):



        #accessing the excel file path from json
        excel_file_path = data["analysis_parameters"]["excel_file_path"]

        #loading the analysis sheet into a dataframe
        analysis_sheet = pd.read_excel(excel_file_path, engine='openpyxl')

        #removing all the rows which has NA values ( basically empty rows after the last sample_id)
        analysis_sheet = analysis_sheet.dropna(how='all')
        return analysis_sheet



#logic for getting the guide sequences(check)
def get_guide_bases(guide):
    my_lims = lims_utils.LimsProcess()
    my_entity = my_lims.get_entity_by_registry_id(guide)
    guide_seq = my_entity.json()['bases']
    return guide_seq


def ProcessAnalysisSheet(data,analysis_sheet):

    samples = []
    #path = analysis_sheet["run_id"][0]
    #done this to get the PRJNA Number or the REQ id from the path so that we dont hardcode it
    #project_id = path.split('/')[-1]
    project_id = data['analysis_parameters']["req_id"]
    bucket_name = data['analysis_parameters']["bucket_name"]




    for _, row in analysis_sheet.iterrows():
        sample_name = row['sample_id']

        guide_1 = row['guide_1']
        guide1_seq = get_guide_bases(guide_1)
        guide_2 = row['guide_2']
        guide2_seq = get_guide_bases(guide_2)

        #creating a dictionery known as samples
        samples.append({
            "sample_name": sample_name,
            "project_id": project_id,
            "guide1": guide1_seq,
            "guide2": guide2_seq,
            #did not give s3 bam path in Excel thought this is much easier
            "s3_bam_path": f"{bucket_name}/{project_id}/{sample_name}/alignment/{sample_name}.target.bam"
        })

    return samples



def main():

    #output directory for the downloaded BAM and generated r1 and r2 fastq files
    output_dir = "output_fastq"
    os.makedirs(output_dir, exist_ok=True)

    config = load_json_file("test.json")
    analysis_sheet = download_analysis_sheet(config)
    samples = ProcessAnalysisSheet(config,analysis_sheet)


   #getting the values from dictionery
    for sample in samples:
        sample_name = sample["sample_name"]
        s3_bam_path = sample["s3_bam_path"]
        local_bam_path = os.path.join(output_dir, f"{sample_name}.bam")

        try:
            # Add check=True to raise an exception if the command fails
            subprocess.run(['aws', 's3', 'cp', s3_bam_path, local_bam_path], check=True)
            print(f"Downloaded {sample_name}")

            # Create the full shell command
            # Fix the path formatting and add proper spacing
            command = f"samtools fastq -1 {os.path.join(output_dir, f'target_{sample_name}_r1.fastq')} \
                          -2 {os.path.join(output_dir, f'target_{sample_name}_r2.fastq')} \
                          -0 /dev/null -s /dev/null -n {local_bam_path}"

            # Execute samtools command
            subprocess.run(command.split(), check=True)
            print(f"Converted {sample_name} to FASTQ")

        except Exception as e:
            print(f"Error processing {sample_name}: {e}")


if __name__ == "__main__":
    main()





