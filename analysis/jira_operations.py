import os
from pathlib import Path
from pipe_utils.jira_utils import JiraActions


jira_ticket = "TEST-16836"
file_folder = "/Users/harshini.muthukumar/Downloads"
file_list = ["REQ4508-001-RNAseq-AnalysisSheet.xlsx"]

jira_obj = JiraActions(jira_issue=jira_ticket)
jira_obj.transition(transition="Start Progress")
for file in file_list:
    file_path = os.path.join(file_folder, file)
    jira_obj.attach_file(attachment=file_path,)
    print(f"Successfully attached {file} to {jira_ticket}")
    jira_obj.transition(transition="Done")
