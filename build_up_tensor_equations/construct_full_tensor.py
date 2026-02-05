
import sympy as sp

from R import R
from ProductTerm import ProductTerm
from build_up_tensor_equations.variations_between_R_and_delta import variations_between_R_and_delta




def construct_all_tensor_terms(order:int):
    # ! without R^-(...) prefactor
    summand_terms = []

    for exponent in range(0, order+1, 2):
        p = r_part(order=order, exponent_of_r=exponent)
        for i in p:
            summand_terms.append(i)

    return summand_terms


def prefactor_expansion(order:int):
    n = sp.symbols('n', integer=True, nonnegative=True)
    term = (-1) ** n * sp.factorial2(2 * n - 1)
    prefactor = term.subs(n, order)
    return prefactor


def r_part(order:int, exponent_of_r:int):
    if exponent_of_r % 2 != 0:
        raise Exception("only even numbers fitting for counter part")
    results = []
    if order < exponent_of_r or order == 0:
        return results

    p_r = variations_between_R_and_delta(amount_of_r=order-exponent_of_r, amount_of_delta=exponent_of_r//2)
    for part in p_r:
        if exponent_of_r > 0:
            r = R("z")
            r.exponent = exponent_of_r
            elements = [r] + part.factors
        else:
            elements = part.factors

        p = ProductTerm()
        p.set_elements(elements=elements, prefactor=prefactor_expansion(order=order-exponent_of_r//2))
        results.append(p)

    return results



#INFO: Herleitung der Funktion r_part:
#
# def r_zero_part(order:int) -> list[ProductTerm]:
#     if order < 1:
#         return []
#     p = variations_between_R_and_delta(amount_of_delta=0, amount_of_r=order)
#     if len(p) != 1:
#         raise Exception("something wrong with function variations_between_R_and_delta")
#     p = p[0]
#     p.prefactor = prefactor_expansion(order=order)
#     return [p]
#
#
# def r_two_part(order:int) -> list[ProductTerm]:
#     results = []
#     if order < 2:
#         return results
#
#     p_r = variations_between_R_and_delta(amount_of_r=order-2, amount_of_delta=1)
#     for part in p_r:
#         r = R("z")
#         r.exponent = 2
#
#         p = ProductTerm()
#         p.set_elements(elements=[r]+part.factors, prefactor=prefactor_expansion(order=order-1))
#         results.append(p)
#
#     return results
#
# def r_four_part(order:int) -> list[ProductTerm]:
#     results = []
#     if order < 4:
#         return results
#
#     p_r = variations_between_R_and_delta(amount_of_r=order-4, amount_of_delta=2)
#     for part in p_r:
#         r = R("z")
#         r.exponent = 4
#
#         p = ProductTerm()
#         p.set_elements(elements=[r]+part.factors, prefactor=prefactor_expansion(order=order-2))
#         results.append(p)
#
#     return results
#
#
# def r_six_part(order:int) -> list[ProductTerm]:
#     results = []
#     if order < 6:
#         return results
#
#     p_r = variations_between_R_and_delta(amount_of_r=order - 6, amount_of_delta=3)
#     for part in p_r:
#         r = R("z")
#         r.exponent = 4
#
#         p = ProductTerm()
#         p.set_elements(elements=[r] + part.factors, prefactor=prefactor_expansion(order=order - 3))
#         results.append(p)
#     return results


