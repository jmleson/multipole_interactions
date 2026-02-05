from sympy import pretty

from R import R
from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms
from settings import global_symbols


def print_full_tensor(order:int):
    tensor_terms = construct_all_tensor_terms(order=order)

    r_prefactor = R("z")
    r_prefactor.exponent = 2 * order - 1


    string  = "T_" + "".join([pretty(i) for i in global_symbols[:order]] ) + " = "
    string += f"{r_prefactor.to_string()} * [ "
    string += " ".join([term.to_string() for term in tensor_terms])
    string += " ]"

    print(string)
    return string



if __name__ == "__main__":
    print_full_tensor(2)