from pipe_utils.jira_utils import JiraActions
import os
import argparse




class rnaseq_Jira(JiraActions):
    def __init__(self, jira_issue):
        super().__init__(jira_issue=jira_issue)

    def attach_results_to_jira(self, files_to_attach):
        """
        Attach aggregated results CSV files to a JIRA ticket.
        """
        # Check if files exist and attach them
        for file_name in files_to_attach:
            file_path = os.path.join(os.getcwd(), file_name)

            if os.path.exists(file_path):
                try:
                    # Attach the file to JIRA ticket
                    self.attach_file(
                        attachment=file_path,
                        jira_issue=self.jira_issue,
                        filename=file_name
                    )
                    print(f"Successfully attached {file_name} to {self.jira_issue}")
                except Exception as e:
                    print(f"Failed to attach {file_name}: {str(e)}")
            else:
                print(f"File not found: {file_path}")


def main():
    parser = argparse.ArgumentParser(description="Attach results to JIRA ticket")
    parser.add_argument("jira_issue", help="JIRA issue key")
    parser.add_argument("files", type=str, help="Files to attach", nargs='+')
    args = parser.parse_args()
    jira_client = rnaseq_Jira(jira_issue=args.jira_issue)
    jira_client.attach_results_to_jira(args.files)


if __name__ == "__main__":
    main()














import os
from pathlib import Path
from pipe_utils.jira_utils import JiraActions



class rnaseq_Jira(JiraActions):
    def __init__(self, jira_issue):
        super().__init__(jira_issue=jira_issue)

    # Transition with a comment
    def state_transition(self, transition):
        self.transition(
            transition=transition
        )

    def attach_results_to_jira(self, file_path):
        """
        Attach aggregated results CSV files to a JIRA ticket.
        """
        # Check if files exist and attach them
        filename = str(Path(file_path).name)

        if os.path.exists(file_path):
            try:
                # Attach the file to JIRA ticket
                self.attach_file(
                    attachment=file_path,
                    jira_issue=self.jira_issue,
                    filename=str(filename)
                )
                print(f"Successfully attached {filename} to {self.jira_issue}")
            except Exception as e:
                print(f"Failed to attach {filename}: {str(e)}")
        else:
            print(f"File not found: {file_path}")


def test():
    jira_ticket = "TEST-16836"
    file_folder = "/Users/harshini.muthukumar/Downloads"
    file_list = ["REQ4508-001-RNAseq-AnalysisSheet.xlsx"]

    jira_obj = rnaseq_Jira(jira_issue=jira_ticket)
    jira_obj.state_transition(transition="In Progress")
    for file in file_list:
        file_path = os.path.join(file_folder, file)
        jira_obj.attach_results_to_jira(file_path=file_path)
        print(f"Successfully attached {file} to {jira_ticket}")
    jira_obj.state_transition(transition="Done")


if __name__ == "__main__":
    test()