
from sympy import pretty
import sympy as sp

from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms, prefactor_expansion
from mathematics_of_terms.add_product_terms import add_product_terms, addable_product_terms
from mathematics_of_terms.sum_up_list import sum_up_list
from settings import global_symbols, greek_names
from term_components.Index import Index
from term_components.Multipole import MultipoleMoment
from term_components.R import R
from term_components.get_product_terms_from_greek_sum import get_product_terms_from_greek_sum_duplicateMultipoles, \
    get_product_terms_from_greek_sum_all


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
        self.prefactor_expansion = 1/( prefactor_expansion(order=multipole_order_1) * prefactor_expansion(order=multipole_order_2) )
        self.set_r_prefactor()

    def get_orders(self):
        return (len(self.multipole1.indices), len(self.multipole2.indices))

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
            string += f"[ "
            string += " ".join([term.to_string() for term in self.tensor_terms if not term.prefactor == 0 ])
            string += " ]"
        else:
            string += "+ 0"
        return string

    def latex_width(self):
        width = 100
        if len(self.tensor_terms) <= 5:
            width = 1
        if (self.get_orders() == (1, 1)
                or self.get_orders() == (1, 3) or self.get_orders() == (3, 1)
                or self.get_orders() == (2, 1) or self.get_orders() == (1, 2)
        ):
            width = 1 # always break since derivation isn't too long
        if self.get_orders() == (1,4) or self.get_orders() == (4,1):
            width = 100# 105 is still too big
        if self.get_orders() == (2,3) or self.get_orders() == (3,2):
            width = 105# 110 is still too big
        if self.get_orders() == (2,4) or self.get_orders() == (4,2):
            width = 110
        return width


    def full_latex_tensor(self) -> str:
        string = self.get_name(latex=True) + " = "
        if len(self.tensor_terms) > 0 and self.prefactor_expansion != 0:
            if self.r_prefactor.exponent != 0:
                string += f"{self.r_prefactor.to_latex()}" + r" \cdot "
            if self.prefactor_expansion != 1:
                string += f"{sp.latex(self.prefactor_expansion)}" + r" \cdot "
            if len(self.tensor_terms) > 0 and len([t for t in self.tensor_terms if t.prefactor != 0]) > 0:
                string += r" \left[ \begin{array}{c} "+ "\n \t\t"
                s = " \n " + r"\\[0.2em] " + "\t\t "
                terms = [term.to_latex() for term in self.tensor_terms if not term.prefactor == 0 ]



                # else:
                line = ""
                for i in range(len(terms)):
                    if len(line) < self.latex_width():
                        line += terms[i]
                    else:
                        string += line + s
                        line = terms[i]
                string += line
                string += "\n" + r" \end{array}\right] "
            else:
                string += r"\left[ + 0 \right]"
        else:
            string += "+ 0"
        return string


    def add_up(self):
        if len(self.tensor_terms) == 0:
            return
        result = sum_up_list(self.tensor_terms, addable_function=addable_product_terms, adding_function=add_product_terms )
        if len(result) > len(self.tensor_terms):
            raise Exception("should never happen")
        self.tensor_terms = result


    def find_unresolved_multipoles(self, function):
        if function != get_product_terms_from_greek_sum_duplicateMultipoles and function != get_product_terms_from_greek_sum_all:
            raise Exception("undefined behavior")
        new_terms = []
        changed_anything = False
        for p in self.tensor_terms:
            result = function(p)
            changed_anything = changed_anything or len(result) > 1
            for term in result:
                new_terms.append(term)
        if changed_anything:
            self.tensor_terms = new_terms
        return changed_anything


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

        self.find_unresolved_multipoles(function=get_product_terms_from_greek_sum_duplicateMultipoles)
        self.find_unresolved_multipoles(function=get_product_terms_from_greek_sum_all)
        self.add_up()
        self.clean_up()



    def add_step(self, latex_doc, description, term_clean_up:bool = True, prior_equation:str="") -> str:
            basic_equation = self.full_latex_tensor()
            if basic_equation == prior_equation:
                return basic_equation

            latex_doc.append(r"\section*{" + description + ":}")

            cleaned_up_equation, added_up_equation = None, None
            if term_clean_up:
                for term in self.tensor_terms:
                    term.clean_up()
                cleaned_up_equation = self.full_latex_tensor()
                self.add_up()
                added_up_equation = self.full_latex_tensor()

            last_equation = basic_equation
            if not term_clean_up or (cleaned_up_equation == basic_equation and added_up_equation == basic_equation):
                latex_doc.append(r"\begin{equation}")
                latex_doc.append(basic_equation)
                latex_doc.append(r"\end{equation}")

            else:
                latex_doc.append(r"\begin{subequations}\begin{gather}")
                latex_doc.append(basic_equation)

                if cleaned_up_equation != basic_equation:
                    latex_doc.append(r"\\ " + r"\xRightarrow{\text{cleaned up}}")
                    latex_doc.append(cleaned_up_equation[cleaned_up_equation.find("=") + 1:])
                    last_equation = cleaned_up_equation
                if added_up_equation != basic_equation and added_up_equation != cleaned_up_equation:
                    latex_doc.append(r"\\ " + r"\xRightarrow{\text{added up}}")
                    latex_doc.append(added_up_equation[added_up_equation.find("=") + 1:])
                    last_equation = added_up_equation

                latex_doc.append(r"\end{gather}\end{subequations}")
            latex_doc.append("\n")
            return last_equation

    def simplify_in_latex_steps(self):
        filename = f"tensor_simplification_{len(self.multipole1.indices)}-{len(self.multipole2.indices)}.tex"
        # Start the LaTeX document
        latex_doc = []
        latex_doc.append(r"\documentclass{article}")
        latex_doc.append(r"\usepackage[top=2cm,left=2cm,bottom=2cm,right=2cm]{geometry}")
        latex_doc.append(r"\usepackage{amsmath}% for equation environment")
        latex_doc.append(r"\usepackage{mathtools}% for xRightarrow")
        latex_doc.append(r"\allowdisplaybreaks % allows breaks in subequations")
        latex_doc.append(r"\begin{document}")



        prior_equation = self.add_step(latex_doc=latex_doc, description=f"Interacting Multipoles of Rank {len(self.multipole1.indices)} and {len(self.multipole2.indices)}",
                 term_clean_up=False)

        for term in self.tensor_terms:
            term.simplify_delta()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Simplify deltas", prior_equation=prior_equation)

        for term in self.tensor_terms:
            term.simplify_R()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Simplify R", prior_equation=prior_equation)

        for max_order, name in [(1, "Dipoles"),
                                (2, "Quadrupoles"),
                                (3, "Oktopoles"),
                                (4, "Hexadecapoles")]:
            for term in self.tensor_terms:
                term.simplify_Multipole(max_order=max_order)
            prior_equation = self.add_step(latex_doc=latex_doc, description=f"Simplify {name}", prior_equation=prior_equation)

        for term in self.tensor_terms:
            term.use_traceless_conditions()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Exploit Traceless Conditions", prior_equation=prior_equation)

        for term in self.tensor_terms:
            term.reduce_indices()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Reduce Indices", prior_equation=prior_equation)

        change = self.find_unresolved_multipoles(function=get_product_terms_from_greek_sum_duplicateMultipoles)
        if change:
            prior_equation = self.add_step(latex_doc=latex_doc, description="Multiply Out Sum (Duplikate Multipoles)", prior_equation=prior_equation)
        change = self.find_unresolved_multipoles(function=get_product_terms_from_greek_sum_all)
        if change:
            prior_equation = self.add_step(latex_doc=latex_doc, description="Multiply Out Sum", prior_equation=prior_equation)

        self.clean_up()
        self.add_step(latex_doc=latex_doc, description="Combining everything", prior_equation=prior_equation)

        latex_doc.append(r"\end{document}")
        with open(filename, "w") as f:
            f.write("\n".join(latex_doc))
        print(f"LaTeX document written to {filename}")
