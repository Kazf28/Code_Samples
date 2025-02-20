import subprocess

task_html_path = "/SPECIFY_PATH/generated_task.html"
solution_html_path = "/SPECIFY_PATH/generated_solution.html"
task_pdf_path = "/SPECIFY_PATH/generated_task.pdf"
solution_pdf_path = "/SPECIFY_PATH/generated_solution.pdf"

subprocess.run(["wkhtmltopdf", "file://" + task_html_path, task_pdf_path])
subprocess.run(["wkhtmltopdf", "file://" + solution_html_path, solution_pdf_path])