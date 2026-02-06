from sympy import pretty
import sympy as sp

from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms, prefactor_expansion
from mathematics_of_terms.add_product_terms import add_product_terms, addable_product_terms
from mathematics_of_terms.sum_up_list import sum_up_list
from settings import global_symbols, greek_names
from term_components.Index import Index
from term_components.Multipole import MultipoleMoment
from term_components.R import R


class MultipoleInteraction:

    def __init__(self, multipole_order_1:int, multipole_order_2:int):

        self.multipole1 = MultipoleMoment(greek_names[:multipole_order_1])
        self.multipole2 = MultipoleMoment(greek_names[multipole_order_1:multipole_order_1+multipole_order_2])

        self.order = len(self.multipole1.indices) + len(self.multipole2.indices)

        self.tensor_terms = []
        for term in construct_all_tensor_terms(order=self.order):
            term.factors.append(  MultipoleMoment(greek_names[:multipole_order_1]) )
            term.factors.append( MultipoleMoment(greek_names[multipole_order_1:multipole_order_1+multipole_order_2]) )
            self.tensor_terms.append(term)

        # initial values (! get changed by clean up):
        self.prefactor_expansion = 1/prefactor_expansion(order=self.order-1)
        self.set_r_prefactor()

    def set_for_testing(self, tensor_terms, order):
        self.multipole1 = None
        self.multipole2 = None
        self.order = order
        self.tensor_terms = tensor_terms

    def set_r_prefactor(self):
        self.r_prefactor = R("z")
        self.r_prefactor.exponent = - (2 * self.order + 1)

    def get_name(self, latex:bool=False):
        if latex:
            string = r"T_{" + "".join([Index(str(i)).to_latex() for i in global_symbols[:self.order]])
            string += r"} " + self.multipole1.to_latex() + " " + self.multipole2.to_latex()
        else:
            string = "T_" + "".join([pretty(i) for i in global_symbols[:self.order]])
            string += " " + self.multipole1.to_string() + " " + self.multipole2.to_string()
        return string

    def full_string_tensor(self):
        string = self.get_name() + " = "
        if len(self.tensor_terms) > 0 and self.prefactor_expansion != 0:
            if self.r_prefactor.exponent != 0:
                string += f"{self.r_prefactor.to_string()} * "
            if self.prefactor_expansion != 1:
                string += f"{self.prefactor_expansion} * "
            string += f" [ "
            string += " ".join([term.to_string() for term in self.tensor_terms if not term.prefactor == 0 ])
            string += " ]"
        else:
            string += "+ 0"
        return string

    def full_latex_tensor(self):
        string = self.get_name(latex=True) + " = "
        if len(self.tensor_terms) > 0 and self.prefactor_expansion != 0:
            if self.r_prefactor.exponent != 0:
                string += f"{self.r_prefactor.to_latex()}" + r" \cdot "
            if self.prefactor_expansion != 1:
                string += f"{sp.latex(self.prefactor_expansion)}" + r" \cdot "
            if len(self.tensor_terms) == 0:
                string += r" \left[ \begin{array}{c} "+ "\n \t\t"
                s = " \n " + r"\\ " + "\t\t "
                string += s.join([term.to_latex() for term in self.tensor_terms if not term.prefactor == 0 ])
                string += "\n" + r" \end{array}\right] "
            else:
                string += r"\left[ + 0 \right]"
        else:
            string += "+ 0"
        return string


    def add_up(self):
        result = sum_up_list(self.tensor_terms, addable_function=addable_product_terms, adding_function=add_product_terms )
        if len(result) > len(self.tensor_terms):
            raise Exception("should never happen")
        self.tensor_terms = result


    def clean_up(self):
        tensor_terms = []
        for term in self.tensor_terms:
            if term.prefactor != 0:
                tensor_terms.append(term)
        self.tensor_terms = tensor_terms

        if self.prefactor_expansion != 1:
            for term in self.tensor_terms:
                term.prefactor *= self.prefactor_expansion
        self.prefactor_expansion = 1

        if self.r_prefactor.exponent != 0:
            for term in self.tensor_terms:
                r = R("z")
                r.exponent = self.r_prefactor.exponent
                term.factors.append(r)
                term.clean_up()
            self.r_prefactor.exponent = 0

    def simplify(self):
        for term in self.tensor_terms:
            term.simplify_delta()
            term.simplify_R()
            term.simplify_Multipole()

            term.use_traceless_conditions()

            term.clean_up()
            term.reduce_indices()
            term.clean_up()

        self.add_up()
        self.clean_up()


    def simplify_in_latex_steps(self):
        filename = f"tensor_simplification_{len(self.multipole1.indices)}-{len(self.multipole2.indices)}.tex"
        # Start the LaTeX document
        latex_doc = []
        latex_doc.append(r"\documentclass{article}")
        latex_doc.append(r"\usepackage{amsmath}% for equation environment")
        latex_doc.append(r"\usepackage{mathtools}% for xRightarrow")
        latex_doc.append(r"\begin{document}")

        # Helper to add an equation with a description
        def add_step(description, term_clean_up:bool=True):
            latex_doc.append(r"\section*{" + description + "}")
            if not term_clean_up:
                latex_doc.append(r"\begin{equation}")
                latex_doc.append(self.full_latex_tensor())
                latex_doc.append(r"\end{equation}")
            else:
                latex_doc.append(r"\begin{subequations}\begin{gather}")
                latex_doc.append(self.full_latex_tensor())
                for term in self.tensor_terms:
                    term.clean_up()
                latex_doc.append(r"\\ " + r"\xRightarrow{\text{cleaned up:}}")
                full_equation = self.full_latex_tensor()
                latex_doc.append(full_equation[full_equation.find("=")+1:])
                latex_doc.append(r"\end{gather}\end{subequations}")
            latex_doc.append("\n")

        add_step(f"Interacting Multipoles of Rank {len(self.multipole1.indices)} and {len(self.multipole2.indices)}:", term_clean_up=False)

        for term in self.tensor_terms:
            term.simplify_delta()
        add_step("Simplify deltas")

        for term in self.tensor_terms:
            term.simplify_R()
        add_step("Simplify R")

        for term in self.tensor_terms:
            term.simplify_Multipole()
        add_step("Simplify Multipoles")

        for term in self.tensor_terms:
            term.use_traceless_conditions()
        add_step("Exploit Traceless Conditions")

        for term in self.tensor_terms:
            term.reduce_indices()
        add_step("Reduce Indices")

        self.add_up()
        add_step("Add up", term_clean_up=False)

        self.clean_up()
        add_step("Combining everything", term_clean_up=False)


        latex_doc.append(r"\end{document}")
        with open(filename, "w") as f:
            f.write("\n".join(latex_doc))
        print(f"LaTeX document written to {filename}")
