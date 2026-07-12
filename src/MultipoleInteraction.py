import sympy as sp

from src.build_up_tensor_equations.construct_full_tensor import prefactor_expansion
from src.latex_helping_functions import save_lines_as_latex_file
from src.mathematics_of_terms.add_product_terms import add_product_terms, addable_product_terms
from src.mathematics_of_terms.sum_up_list import sum_up_list
from src.settings import greek_names
from src.term_components.Multipole import MultipoleMoment
from src.term_components.ProductTerm import ProductTerm
from src.term_components.debug_replacements import get_replacements_according_to_traceless_conditions
from src.term_components.get_product_terms_from_greek_sum import get_product_terms_from_greek_sum_duplicateMultipoles, \
    get_product_terms_from_greek_sum_all
from src.Tensor import Tensor


class MultipoleInteraction:

    def __init__(self, multipole_order_1:int, multipole_order_2:int):

        self.multipole1 = MultipoleMoment(greek_names[:multipole_order_1])
        self.multipole1.molecule = "B"
        self.multipole2 = MultipoleMoment(greek_names[multipole_order_1:multipole_order_1+multipole_order_2])
        self.multipole2.molecule = "A"

        self.order = len(self.multipole1.indices) + len(self.multipole2.indices)

        self.tensor = Tensor(self.order)
        # add multipoles to tensor terms:
        for t in range(len(self.tensor.terms)):
            m1 = MultipoleMoment(greek_names[:multipole_order_1])
            m1.molecule = "B"
            self.tensor.terms[t].factors.append( m1 )
            m2 = MultipoleMoment(greek_names[multipole_order_1:multipole_order_1+multipole_order_2])
            m2.molecule = "A"
            self.tensor.terms[t].factors.append( m2 )

        # initial values (! get changed by clean up):
        self.prefactor_of_expansion = abs(1 / prefactor_expansion(order=multipole_order_1)) * (1/prefactor_expansion(multipole_order_2))
        self.initial_prefactor = self.prefactor_of_expansion

        self.full_results = []

    def get_orders(self):
        return (len(self.multipole1.indices), len(self.multipole2.indices))

    def get_name(self, latex:bool=False):
        if latex:
            string = "+ " if self.initial_prefactor >= 0 else "- "
            string += sp.latex(abs(self.initial_prefactor)) + r" \cdot "
            string += self.tensor.get_name(latex=latex)
            string += r" " + self.multipole1.to_latex() + " " + self.multipole2.to_latex()
        else:
            string = "+ " if self.initial_prefactor >= 0 else "- "
            string += str(abs(self.initial_prefactor)) + r" * "
            string += self.tensor.get_name(latex=latex)
            string += " " + self.multipole1.to_string() + " " + self.multipole2.to_string()
        return string

    def full_string_tensor(self):
        string = self.get_name() + " = "
        if len(self.tensor.terms) > 0 and self.prefactor_of_expansion != 0:
            if self.tensor.r_prefactor.exponent != 0:
                string += f"{self.tensor.r_prefactor.to_string()} * "
            if self.prefactor_of_expansion != 1:
                string += f"{self.prefactor_of_expansion} * "
            string += f"[ "
            string += " ".join([term.to_string() for term in self.tensor.terms if not term.prefactor == 0 ])
            string += " ]"
        else:
            string += "+ 0"
        return string

    def latex_width(self):
        width = 100
        if len(self.tensor.terms) <= 5:
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
        if len(self.tensor.terms) > 0 and self.prefactor_of_expansion != 0:
            if self.tensor.r_prefactor.exponent != 0:
                string += f"{self.tensor.r_prefactor.to_latex()}" + r" \cdot "
            if self.prefactor_of_expansion != 1:
                string += f"{sp.latex(self.prefactor_of_expansion)}" + r" \cdot "
            if len(self.tensor.terms) > 0 and len([t for t in self.tensor.terms if t.prefactor != 0]) > 0:
                string += r" \left[ \begin{array}{c} "+ "\n \t\t"
                s = " \n " + r"\\[0.2em] " + "\t\t "
                terms = [term.to_latex() for term in self.tensor.terms if not term.prefactor == 0 ]

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
        if len(self.tensor.terms) == 0:
            return
        result = sum_up_list(self.tensor.terms, addable_function=addable_product_terms, adding_function=add_product_terms )
        if len(result) > len(self.tensor.terms):
            raise Exception("should never happen")
        self.tensor.terms = result


    def find_unresolved_multipoles(self, function):
        if function != get_product_terms_from_greek_sum_duplicateMultipoles and function != get_product_terms_from_greek_sum_all:
            raise Exception("undefined behavior")
        new_terms = []
        changed_anything = False
        for p in self.tensor.terms:
            result = function(p)
            changed_anything = changed_anything or len(result) > 1
            for term in result:
                new_terms.append(term)
        if changed_anything:
            self.tensor.terms = new_terms
        return changed_anything


    def clean_up(self):
        self.tensor.clean_up()

        if self.prefactor_of_expansion != 1:
            for term in self.tensor.terms:
                term.prefactor *= self.prefactor_of_expansion
        self.prefactor_of_expansion = 1



    def add_step(self, latex_doc, description, term_clean_up:bool = True, prior_equation:str="", consider_index_sorting:bool=False, importable_tex:bool=False) -> str:
            basic_equation = self.full_latex_tensor()
            if basic_equation == prior_equation:
                return basic_equation


            if prior_equation is None:
                if importable_tex:
                    latex_doc.append(r"\subsubsection{" + description + "}")
                else:
                    latex_doc.append(r"\subsubsection*{" + description + "}")
            else:
                latex_doc.append(r"\textbf{" + description + ":}")

            cleaned_up_equation, added_up_equation, reduced_indices_equation = None, None, None
            if term_clean_up:
                for term in self.tensor.terms:
                    term.clean_up()
                cleaned_up_equation = self.full_latex_tensor()

                for term in self.tensor.terms:
                    term.reduce_indices(consider_index_sorting=consider_index_sorting)
                reduced_indices_equation = self.full_latex_tensor()

                self.add_up()
                added_up_equation = self.full_latex_tensor()

            last_equation = basic_equation
            if not term_clean_up or (cleaned_up_equation == basic_equation and added_up_equation == basic_equation and reduced_indices_equation == basic_equation):
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
                if reduced_indices_equation != basic_equation and reduced_indices_equation != cleaned_up_equation:
                    latex_doc.append(r"\\ " + r"\xRightarrow{\text{reduced indices}}")
                    latex_doc.append(reduced_indices_equation[reduced_indices_equation.find("=") + 1:])
                    last_equation = reduced_indices_equation
                if added_up_equation != basic_equation and added_up_equation != reduced_indices_equation:
                    latex_doc.append(r"\\ " + r"\xRightarrow{\text{added up}}")
                    latex_doc.append(added_up_equation[added_up_equation.find("=") + 1:])
                    last_equation = added_up_equation

                latex_doc.append(r"\end{gather}\end{subequations}")
            latex_doc.append("\n")
            return last_equation




    def simplify_in_latex_steps(self, importable_tex:bool=False , debug:bool=False):
        filename = f"tensor_simplification_{len(self.multipole1.indices)}-{len(self.multipole2.indices)}.tex"
        latex_doc = []

        prior_equation = self.add_step(latex_doc=latex_doc, prior_equation=None, importable_tex=importable_tex,
                                       description=f"Interacting Multipoles of Rank {len(self.multipole1.indices)} and {len(self.multipole2.indices)}",
                 term_clean_up=False)

        for term in self.tensor.terms:
            term.simplify_delta()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Simplify deltas", importable_tex=importable_tex, prior_equation=prior_equation, term_clean_up=True)

        for term in self.tensor.terms:
            term.simplify_R()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Simplify R", importable_tex=importable_tex, prior_equation=prior_equation, term_clean_up=True)


        for max_order, name in [(1, "Dipoles"),
                                (2, "Quadrupoles"),
                                (3, "Oktopoles"),
                                (4, "Hexadecapoles")]:
            for term in self.tensor.terms:
                term.simplify_Multipole(max_order=max_order)
            prior_equation = self.add_step(latex_doc=latex_doc, description=f"Simplify {name}",
                                           importable_tex=importable_tex, prior_equation=prior_equation,
                                           term_clean_up=True, consider_index_sorting=True)

        for term in self.tensor.terms:
            term.use_traceless_conditions()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Exploit Traceless Conditions",
                                       importable_tex=importable_tex, term_clean_up=True,
                                       prior_equation=prior_equation, consider_index_sorting=True)


        change = self.find_unresolved_multipoles(function=get_product_terms_from_greek_sum_all)
        if change:
            prior_equation = self.add_step(latex_doc=latex_doc, description="Multiply Out Sum",
                                           importable_tex=importable_tex, prior_equation=prior_equation,
                                           term_clean_up=True, consider_index_sorting=True)

        self.clean_up()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Combining everything",
                                       importable_tex=importable_tex, term_clean_up=False,
                                       prior_equation=prior_equation, consider_index_sorting=True)
        self.full_results.append(prior_equation)
        if debug:
            self.debug_substitute_multipoles_according_to_traceless_conditions("A")
            prior_equation = self.add_step(latex_doc=latex_doc, description="A", term_clean_up=True, importable_tex=importable_tex,
                                           prior_equation=prior_equation, consider_index_sorting=True)
            self.debug_substitute_multipoles_according_to_traceless_conditions("B")
            prior_equation = self.add_step(latex_doc=latex_doc, description="B", term_clean_up=True, importable_tex=importable_tex,
                                           prior_equation=prior_equation, consider_index_sorting=True)
            self.clean_up()
            prior_equation = self.add_step(latex_doc=latex_doc, description="Combining everything",
                                           importable_tex=importable_tex, term_clean_up=False,
                                           prior_equation=prior_equation, consider_index_sorting=True)

            self.full_results.append(prior_equation)

        save_lines_as_latex_file(filename=filename, importable_tex=importable_tex, content = latex_doc)

    def simplify_long_way(self, debug=True):
        filename = f"debug_tensor_simplification_{len(self.multipole1.indices)}-{len(self.multipole2.indices)}.tex"
        latex_doc = []

        prior_equation = self.add_step(latex_doc=latex_doc, prior_equation=None,
                                       description=f"Interacting Multipoles of Rank {len(self.multipole1.indices)} and {len(self.multipole2.indices)}",
                                       term_clean_up=False)

        for term in self.tensor.terms:
            term.simplify_delta()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Simplify deltas",
                                       prior_equation=prior_equation, term_clean_up=True)

        for term in self.tensor.terms:
            term.simplify_R()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Simplify R", prior_equation=prior_equation,
                                       term_clean_up=True)

        change = self.find_unresolved_multipoles(function=get_product_terms_from_greek_sum_all)
        if change:
            prior_equation = self.add_step(latex_doc=latex_doc, description="Multiply Out Sum",
                                           prior_equation=prior_equation,
                                           term_clean_up=True, consider_index_sorting=True)

        self.clean_up()
        prior_equation = self.add_step(latex_doc=latex_doc, description="Combining everything", term_clean_up=False,
                                       prior_equation=prior_equation, consider_index_sorting=True)
        self.full_results.append(prior_equation)

        if debug:
            self.debug_substitute_multipoles_according_to_traceless_conditions("A")
            prior_equation = self.add_step(latex_doc=latex_doc, description="A", term_clean_up=True,
                                           prior_equation=prior_equation, consider_index_sorting=True)
            self.debug_substitute_multipoles_according_to_traceless_conditions("B")
            prior_equation = self.add_step(latex_doc=latex_doc, description="B", term_clean_up=True,
                                           prior_equation=prior_equation, consider_index_sorting=True)
            self.clean_up()
            prior_equation = self.add_step(latex_doc=latex_doc, description="Combining everything", term_clean_up=False,
                                           prior_equation=prior_equation, consider_index_sorting=True)

            self.full_results.append(prior_equation)

        save_lines_as_latex_file(filename=filename, importable_tex=False, content=latex_doc)



    def debug_substitute_multipoles_according_to_traceless_conditions(self, molecule:str):
        factors = []
        for term in self.tensor.terms:
            something_replaced = False
            for i in get_replacements_according_to_traceless_conditions(molecule=molecule):
                if i["to_be_replaced"] in term.factors:
                    something_replaced = True
                    basic_elements = [x for x in term.factors if x != i["to_be_replaced"]]
                    for j in i["replacements"]:
                        p = ProductTerm()
                        p.set_elements(basic_elements+[j], prefactor=-term.prefactor)
                        p.sort_elements()
                        factors.append(p)
            if not something_replaced:
                p = ProductTerm()
                elements = term.factors
                p.set_elements(elements, prefactor = term.prefactor)
                factors.append(p)
        self.tensor.terms = factors


    def __eq__(self, other):
        if not isinstance(other, MultipoleInteraction):
            return NotImplemented
        if self.full_latex_tensor() == other.full_latex_tensor():
            return True
        if set(self.full_results) & set(other.full_results):
            return True
        return False
