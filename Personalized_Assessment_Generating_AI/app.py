from flask import Flask, render_template, request
import Main_LLM_File  # Import your script with the AI functions
import subprocess

app = Flask(__name__)


def write_tips():
    with open("teacher_tips.txt") as f:
        file = f.readlines()
        tips_dict = {}
        for line in file:
            col_index = line.find(":")
            concept = line[:col_index].upper()
            tip = line[col_index + 1:].rstrip()
            tips_dict[concept] = tip
        return tips_dict

def save_task_to_html(result_task, result_solution, concept, name):
    
    task_html_path = f"/SPECIFY_PATH/{name}_{concept}_task.html"
    solution_html_path = f"/SPECIFY_PATH{name}_{concept}_solution.html"

    with open(task_html_path, "w") as file:
        file.write(result_task)
    with open(solution_html_path, "w") as file:
        file.write(result_solution)

    # Use Puppeteer to generate PDFs
    task_pdf_path = f"/SPECIFY_PATH{name}_{concept}_task.pdf"
    solution_pdf_path = f"/SPECIFY_PATH{name}_{concept}_solution.pdf"

    subprocess.run(["node", "generate_pdf.js", task_html_path, task_pdf_path])
    subprocess.run(["node", "generate_pdf.js", solution_html_path, solution_pdf_path])



@app.route('/', methods=['GET', 'POST'])
def index():
    result_task, result_solution, task_name, solution_name = "", "", "", ""

    if request.method == 'POST':
        grade = request.form.get('grade')
        concept = request.form.get('concept')
        interest = request.form.get('interest')
        name = request.form.get('name')
        additional = request.form.get('additional')
        task_name = name + "_" + concept + "_task.pdf"
        solution_name = name + "_" + concept + "_solution.pdf"


        # Generate task and solution from Main_LLM_File
        result_task, result_solution = Main_LLM_File.AIgenerate(grade, concept, interest, additional)

        # Save the generated HTML and PDF files
        save_task_to_html(result_task, result_solution, concept, name)

    return render_template('Digital_Prototype.html', result_task=result_task, result_solution=result_solution, task_name = task_name, solution_name = solution_name)


if __name__ == '__main__':
    app.run(debug=True)
