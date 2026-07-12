from src.term_components.ProductTerm import ProductTerm


def add_product_terms(p1:ProductTerm, p2:ProductTerm):
    if p1.prefactor == 0:
        return [p2]
    if p2.prefactor == 0:
        return [p1]
    if not p1.addable(p2):
        return [p1, p2]

    new_p = ProductTerm()
    new_p.set_elements(elements=p1.factors, prefactor=p1.prefactor + p2.prefactor)

    return [new_p]


def addable_product_terms(p1:ProductTerm, p2:ProductTerm):
    return p1.addable(p2)
