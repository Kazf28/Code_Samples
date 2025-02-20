from langchain_openai import OpenAIEmbeddings, ChatOpenAI
from langchain.prompts import ChatPromptTemplate, HumanMessagePromptTemplate, SystemMessagePromptTemplate,AIMessagePromptTemplate
import os
os.environ["OPENAI_API_KEY"] = "SPECIFY_API_KEY"

def AIgenerate(grade, concept, interest, additional):
    llm_task = ChatOpenAI(model_name="gpt-4o-mini", temperature=0.6)

    # Task generation
    task_system_template = (
        "You are a professional learning designer tasked with creating homework for high school mathematics students. The school has hired you to develop assignments that are 1) tailored to students’ interests and 2) apply specific mathematical concepts to real-life situations. \n We will provide you with two inputs: 1) the mathematical concepts to be covered in the homework, and 2) the students' interests. Based on these inputs, your job is to design a homework worksheet that engages struggling students and motivates them to complete the tasks, making it a story-telling format so that it is intuitive for students to grasp the task, without too much text. Make sure the difficulty matches the grade level of the students. Please ensure that the assignment follows a constructionist approach and is incremental in difficulty, avoiding simple calculation problems but also not to go beyond the high school level mathematics. Also, make sure you avoid using mathematical notations like P(x) for the price function, or avoid ambiguous variables like “a” or “t” and instead, replace them with the word of concept it is actually representing, but make sure you include the function in the question. Do not include instructions on submission deadlines. Produce output ONLY as a html code without unnecessary comment outside the <!DOCTYPE html>. Also, add <script src='htps://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script> and make ALL EQUATIONS in latex format. </html>."
    )
    task_system_message_prompt = SystemMessagePromptTemplate.from_template(task_system_template)
    task_human_template = "The inputs are: \n {GRADE} \n {CONCEPT} \n {INTEREST}. \n Additional Information: {ADDITIONAL}"
    task_human_message_prompt = HumanMessagePromptTemplate.from_template(task_human_template)
    task_prompt = ChatPromptTemplate.from_messages([task_system_message_prompt, task_human_message_prompt])
    task_prompt = task_prompt.format_prompt(GRADE=grade, CONCEPT=concept, INTEREST=interest, ADDITIONAL = additional).to_messages()
    result_task = llm_task(task_prompt).content

    # Solution generation
    llm_solution = ChatOpenAI(model_name="gpt-4o-mini", temperature=0.6)
    solution_system_template = "You are a professional mathematics teacher that gives solutions for inputted assignments. Bullet point key steps that has to be taken when students solve the problem. Produce output ONLY as a html code without unnecessary comment outside the <!DOCTYPE html>. Also, add <script src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script> and make ALL EQUATIONS in latex format. </html>."
    solution_system_message_prompt = SystemMessagePromptTemplate.from_template(solution_system_template)
    solution_human_template = "The input assignment is: {TASK}"
    solution_human_message_prompt = HumanMessagePromptTemplate.from_template(solution_human_template)
    solution_prompt = ChatPromptTemplate.from_messages([solution_system_message_prompt, solution_human_message_prompt])
    solution_prompt = solution_prompt.format_prompt(TASK=result_task).to_messages()
    result_solution = llm_solution(solution_prompt).content

    return result_task, result_solution

