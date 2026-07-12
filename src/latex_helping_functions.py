

def get_preamble():
    latex_doc = []
    latex_doc.append(r"\documentclass{article}")
    latex_doc.append(r"\usepackage[top=2cm,left=2cm,bottom=2cm,right=2cm]{geometry}")
    latex_doc.append(r"\usepackage{amsmath}% for equation environment")
    latex_doc.append(r"\usepackage{mathtools}% for xRightarrow")
    latex_doc.append(r"\allowdisplaybreaks % allows breaks in subequations")
    latex_doc.append(r"\usepackage{breqn}% dmath environment (needed for very long equations)")
    latex_doc.append(r"\begin{document}")

    return latex_doc


def save_lines_as_latex_file(filename:str, content:list[str], importable_tex:bool=False):
    latex_doc = [] if importable_tex else get_preamble()

    latex_doc += content

    if not importable_tex:
        latex_doc.append(r"\end{document}")
    with open(filename, "w") as f:
        f.write("\n".join(latex_doc))
    print(f"LaTeX document written to {filename}", flush=True)

