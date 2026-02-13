from sympy import pretty

from term_components.R import R
from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms
from settings import global_symbols


def full_string_tensor(order:int):
    tensor.terms = construct_all_tensor_terms(order=order)

    r_prefactor = R("z")
    r_prefactor.exponent = - (2 * order - 1)

    string  = "T_" + "".join([pretty(i) for i in global_symbols[:order]] ) + " = "
    string += f"{r_prefactor.to_string()} * [ "
    string += " ".join([term.to_string() for term in tensor_terms])
    string += " ]"

    return string



def full_latex_tensor(order:int):
    tensor.terms = construct_all_tensor_terms(order=order)

    r_prefactor = R("z")
    r_prefactor.exponent = - (2 * order - 1)

    string  = "T_{" + "".join([r"\ "[0]+str(i) for i in global_symbols[:order]] ) + "} = "
    string += f"{r_prefactor.to_latex()} " + r"\cdot \left[ "
    string += " ".join([term.to_latex() for term in tensor_terms])
    string += r" \right]"

    return string


if __name__ == "__main__":
    print(full_string_tensor(2))

    print( full_latex_tensor(2) )