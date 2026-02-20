
from sympy import pretty

from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms
from latex_helping_functions import save_lines_as_latex_file
from settings import global_symbols, greek_names
from term_components.Index import Index
from term_components.R import R


class Tensor:

    def __init__(self, order:int):
        assert order <= len(greek_names)

        self.order = order
        self.terms = construct_all_tensor_terms(order)

        self.latex_width = 200
        self.latex_max_number_in_array = 50

        self.r_prefactor = R("z")
        self.r_prefactor.exponent = - (2 * self.order + 1)


    def get_name(self, latex: bool = False):
        string = "T"
        if self.order > 0:
            if latex:
                string += r"_{" + "".join([Index(str(i)).to_latex() for i in global_symbols[:self.order]]) + r"}"
            else:
                string += "_" + "".join([pretty(i) for i in global_symbols[:self.order]])
        return string


    def full_string_tensor_basics(self):
        string = self.get_name() + " = "
        if self.r_prefactor.exponent != 0:
            string += f"{self.r_prefactor.to_string()}"
        if len(self.terms) > 0 :
            string += f" * [ "
            string += " ".join([term.to_string() for term in self.terms if not term.prefactor == 0])
            string += " ]"
        return string.replace("R_z", "R")

    def clean_up(self):
        tensor_terms = []
        for term in self.terms:
            if term.prefactor != 0:
                tensor_terms.append(term)
        self.terms = tensor_terms

        if self.r_prefactor.exponent != 0:
            for term in self.terms:
                r = R("z")
                r.exponent = self.r_prefactor.exponent
                term.factors.append(r)
                term.clean_up()
            self.r_prefactor.exponent = 0

    # def latex_dmath(self) -> str:
    #     """
    #     Erzeugt LaTeX-Code für eine lange Tensor-Gleichung über mehrere dmath-Blöcke
    #     mit offener Klammer am Anfang und geschlossener Klammer am Ende.
    #     Ideal für sehr viele Terme, ggf. über Seitenumbrüche.
    #     """
    #     # Präfix der Gleichung: Tensorname und evtl. R-Vorfaktor
    #     string = self.get_name(latex=True) + " = "
    #     if self.r_prefactor.exponent != 0:
    #         string += f"{self.r_prefactor.to_latex()} "
    #
    #     # Alle Terme, die nicht null sind
    #     terms = [term.to_latex() for term in self.terms if term.prefactor != 0]
    #
    #     if not terms:
    #         string += "+ 0"
    #         return string.replace("R_{z}", "R")
    #
    #     # Blockgröße definieren (maximale Länge pro dmath-Block)
    #     max_block_len = 5  # Anzahl Terme pro dmath-Block, passt man je nach Länge an
    #
    #     # Beginne subequations
    #     latex_blocks = [r"\begin{subequations}"]
    #
    #     # Offene Klammer am ersten Block
    #     first_block_terms = terms[:max_block_len]
    #     first_block = r"\begin{dmath}" + "\n" + string + r"\Big[" + "\n"
    #     first_block += "\n".join(first_block_terms) + "\n" + r"\end{dmath}"
    #     latex_blocks.append(first_block)
    #
    #     # Mittlere Blöcke (dmath* ohne Nummer)
    #     for i in range(max_block_len, len(terms), max_block_len):
    #         block_terms = terms[i:i + max_block_len]
    #         block = r"\begin{dmath*}" + "\n" + " + ".join(block_terms) + "\n" + r"\end{dmath*}"
    #         latex_blocks.append(block)
    #
    #     # Schließende Klammer am letzten Block
    #     latex_blocks[-1] = latex_blocks[-1][:-len(r"\end{dmath*")] + "\ "[0] + r"Big]"+r"\end{dmath*}"
    #
    #     # Ende subequations
    #     latex_blocks.append(r"\end{subequations}")
    #
    #     return "\n\n".join(latex_blocks).replace("R_{z}", "R").replace(r"\\Big]", "\Big]")

    def full_latex_tensor_basics(self) -> str:
        string = self.get_name(latex=True) + " = "
        if self.r_prefactor.exponent != 0:
            string += f"{self.r_prefactor.to_latex()}"
            if len(self.terms) > 0 and len([t for t in self.terms if t.prefactor != 0]) > 0:
                string += r" \cdot "
                if len(self.terms) < self.latex_max_number_in_array:
                    string += r" \left[ \begin{array}{c} " + "\n \t\t"
                else:
                    string += "["+" \n " + r"\\[0.2em] \notag "
                s = " \n " + r"\\[0.2em] \notag " + "\t\t "
                terms = [term.to_latex() for term in self.terms if not term.prefactor == 0 ]

                # else:
                line = ""
                for i in range(len(terms)):
                    if len(line) < self.latex_width:
                        line += terms[i]
                    else:
                        string += line + s
                        line = terms[i]
                string += line
                if len(self.terms) < self.latex_max_number_in_array:
                    string += "\n" + r" \end{array}\right] "
                else:
                    string += r"\\[0.2em] \notag " + "]"
        else:
            string += "+ 0"
        return string.replace("R_{z}", "R")

    def to_tex(self, importable_tex:bool=False):
        filename = f"tensor_definition_order{self.order}"
        import_line = r"\input{" + filename + r"}"
        latex_doc = []

        latex_doc.append(r"\begin{subequations}\begin{gather}")
        latex_doc.append( f"{self.full_latex_tensor_basics()}" )
        latex_doc.append(r"\end{gather}\end{subequations}")
        save_lines_as_latex_file(filename=filename+".tex", importable_tex=importable_tex, content=latex_doc)
        return import_line



if __name__ == '__main__':

    t = Tensor(order=6)
    print(t.full_string_tensor_basics())

    t.to_tex()
