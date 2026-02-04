
import sympy as sp

from Delta import Delta
from R import R
from ProductTerm import ProductTerm
from settings import global_symbols


def prefactor_expansion(order:int):
    n = sp.symbols('n', integer=True, nonnegative=True)
    term = (-1) ** n * sp.factorial2(2 * n - 1)
    prefactor = term.subs(n, order)
    return prefactor



def r_zero_part(order:int):
    terms = []
    for i in range(order):
        symbol = global_symbols[i]
        r = R(index_sign = str(symbol))
        terms.append(r)

    p = ProductTerm()
    p.set_elements(elements = terms, prefactor = prefactor_expansion(order=order))

    return p


def r_two_part(order:int):
    if order <= 1:
        return ProductTerm()
    r = R("z")
    r.exponent = 2

    terms = [r]
    i = 0
    while i < order -2 :
        symbol = global_symbols[i]
        r = R(index_sign = str(symbol))
        terms.append(r)
        i += 1
    d = Delta(str(global_symbols[order - 1]), str(global_symbols[order - 2]))
    terms.append(d)

    p = ProductTerm()
    p.set_elements(elements=terms, prefactor=prefactor_expansion(order=order-1))

    return p




def construct_full_tensor(order:int):
    summand_terms = []
    if order >= 1:
        p = r_zero_part(order)
        summand_terms.append(p)

    if order >= 2:
        p = r_one_part(order)
        summand_terms.append(p)

    return summand_terms