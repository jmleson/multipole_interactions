from sympy import pretty

from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms, prefactor_expansion
from mathematics_of_terms.add_product_terms import add_product_terms, addable_product_terms
from mathematics_of_terms.sum_up_list import sum_up_list
from settings import global_symbols, greek_names
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

    def get_name(self):
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
            term.clean_up()

        self.add_up()
        self.clean_up()



