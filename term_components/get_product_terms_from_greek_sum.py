import copy

from term_components.Delta import Delta
from term_components.Index import Index
from term_components.Multipole import MultipoleMoment
from term_components.ProductTerm import ProductTerm
from term_components.R import R



def get_product_terms_from_greek_sum_all(p:ProductTerm) -> list[ProductTerm]:
    for m in p.factors:
        if isinstance(m, MultipoleMoment):
            greek_indices = m.get_greek_indices()
            if len(greek_indices) > 0:
                # Expand ONE level
                expanded_terms = multiply_out_greek_indices(p, greek_indices)

                # Recurse on each expanded term
                results = []
                for term in expanded_terms:
                    results.extend(get_product_terms_from_greek_sum_all(term))

                return results
        if isinstance(m, Delta):
            raise Exception("delta is supposed to be cleaned beforehand")
        if isinstance(m, R) and not m.index.is_coordinate():
            raise Exception("R is supposed to be set to coordinate beforehand")
    return [p]



def get_product_terms_from_greek_sum_duplicateMultipoles(p:ProductTerm) -> list[ProductTerm]:
    duplicates_included, greek_indices = p.includes_duplicates()
    if not duplicates_included:
        return [p]
    return multiply_out_greek_indices(p, greek_indices=greek_indices)



def multiply_out_greek_indices(p:ProductTerm, greek_indices:list[Index]) -> list[ProductTerm]:
    if len(greek_indices) == 0:
        return [p]

    new_terms = []
    for index in greek_indices:
        # replace by sum of x, y and z variants
        for xyz in ["x", "y", "z"]:
            p_new = ProductTerm()
            to_be_replaced, replacement = str(index.index), xyz

            factors = []
            for i in p.factors:
                if isinstance(i, Delta):
                    index1_str = str(i.index1.index) if str(i.index1.index) != to_be_replaced else xyz
                    index2_str = str(i.index2.index) if str(i.index2.index) != to_be_replaced else xyz
                    new_term = Delta(index1_str, index2_str)
                elif isinstance(i, R):
                    index_sign = str(i.index.index) if str(i.index.index) != to_be_replaced else xyz
                    new_term = R(index_sign)
                    new_term.exponent = i.exponent
                elif isinstance(i, MultipoleMoment):
                    indices = [str(idx.index) if str(idx.index) != to_be_replaced else xyz for idx in i.indices]
                    # print("set by indices", indices)
                    new_term = MultipoleMoment(indices=indices)
                    new_term.molecule = i.molecule
                    # print("->", new_term.to_string())
                else:
                    raise Exception("unknown type of term in ProductTerm")
                factors.append(new_term)
            p_new.set_elements(factors, prefactor=copy.deepcopy(p.prefactor))
            p_new.sort_elements()
            # print("p", p_new.to_string())

            new_terms.append(p_new)
        break # only resolve 1st index

    if len(greek_indices) > 1:
        result = []
        for term in new_terms:
            result.extend(get_product_terms_from_greek_sum_duplicateMultipoles(term))
        return result
    return new_terms