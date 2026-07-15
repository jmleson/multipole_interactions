import sympy as sp


from src.mathematics_of_terms.add_product_terms import addable_product_terms, add_product_terms
from src.mathematics_of_terms.sum_up_list import sum_up_list
from src.term_components.Multipole import MultipoleMoment
from src.term_components.ProductTerm import ProductTerm
from src.term_components.debug_replacements import get_replacements_according_to_traceless_conditions


def debug_substitute_multipoles_according_to_traceless_conditions(terms:list, molecule: str):
    factors = []
    for term in terms:
        something_replaced = False
        for i in get_replacements_according_to_traceless_conditions(molecule=molecule):
            if i["to_be_replaced"] in term.factors:
                something_replaced = True
                basic_elements = [x for x in term.factors if x != i["to_be_replaced"]]
                for j in i["replacements"]:
                    p = ProductTerm()
                    p.set_elements(basic_elements + [j], prefactor=-term.prefactor)
                    p.sort_elements()
                    factors.append(p)
        if not something_replaced:
            p = ProductTerm()
            elements = term.factors
            p.set_elements(elements, prefactor=term.prefactor)
            factors.append(p)
    terms = factors
    return terms



def help_print(result:list):
    for p in result:
        print(p.to_string(), end="\t")
    print()



### spherical result:
sp_term_1 = ProductTerm()
m1 = MultipoleMoment(["y", "z", "z"])
m1.molecule = "A"
m2 = MultipoleMoment(["y", "z", "z"])
m2.molecule = "B"
sp_term_1.set_elements([m1, m2], prefactor = sp.sympify("+15 * 3/2"))


m3 = MultipoleMoment(["x", "x", "y"])
m3.molecule = "A"
m4 = MultipoleMoment(["y", "y", "y"])
m4.molecule = "A"
m5 = MultipoleMoment(["x", "x", "y"])
m5.molecule = "B"
m6 = MultipoleMoment(["y", "y", "y"])
m6.molecule = "B"

sp_term_2 = ProductTerm()
sp_term_2.set_elements([m3, m5], prefactor = sp.sympify("+9/10"))

sp_term_3 = ProductTerm()
sp_term_3.set_elements([m3, m6], prefactor = sp.sympify("-3/10"))

sp_term_4 = ProductTerm()
sp_term_4.set_elements([m4, m5], prefactor = sp.sympify("-3/10"))

sp_term_5 = ProductTerm()
sp_term_5.set_elements([m4, m6], prefactor = sp.sympify("+1/10"))

result_1 = [sp_term_1, sp_term_2, sp_term_3, sp_term_4, sp_term_5]

### cartesian result:
c_term_1 = ProductTerm()
c_term_2 = ProductTerm()
c_term_3 = ProductTerm()

m1 = MultipoleMoment(["y", "z", "z"])
m1.molecule = "A"
m2 = MultipoleMoment(["y", "z", "z"])
m2.molecule = "B"

m3 = MultipoleMoment(["x", "x", "y"])
m3.molecule = "A"
m4 = MultipoleMoment(["x", "x", "y"])
m4.molecule = "B"

m5 = MultipoleMoment(["y", "y", "y"])
m5.molecule = "A"
m6 = MultipoleMoment(["y", "y", "y"])
m6.molecule = "B"

c_term_1.set_elements([m2, m1], prefactor = sp.sympify("+111/5"))
c_term_2.set_elements([m4, m3], prefactor = sp.sympify("+6/5"))
c_term_3.set_elements([m6, m5], prefactor = sp.sympify("+2/5"))#TODO 6/5 statt 2/5

results_2 = [c_term_1, c_term_2, c_term_3]






#######################################





def run(result):
    help_print(result)
    result = sum_up_list(result, addable_function=addable_product_terms, adding_function=add_product_terms)
    # help_print(result)

    result = debug_substitute_multipoles_according_to_traceless_conditions(terms=result, molecule="A")
    # help_print(result)
    result = sum_up_list(result, addable_function=addable_product_terms, adding_function=add_product_terms)
    # help_print(result)
    result = debug_substitute_multipoles_according_to_traceless_conditions(terms=result, molecule="B")
    # help_print(result)
    result = sum_up_list(result, addable_function=addable_product_terms, adding_function=add_product_terms)
    help_print(result)

print("spherical:")
run(result_1)
print()
print("cartesian:")
run(results_2)

