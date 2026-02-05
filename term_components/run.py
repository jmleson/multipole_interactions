from term_components.Delta import Delta
from term_components.Multipole import MultipoleMoment
from term_components.ProductTerm import ProductTerm
from term_components.R import R


#### CASE 1 = Dipole-Dipole-Interaction ####
p1 = ProductTerm()
elements = [
    R("alpha"),
    R("beta"),
    MultipoleMoment(["alpha"]),
    MultipoleMoment(["beta"])
]
p1.set_elements(elements, prefactor=3)
p1.simplify_delta()
p1.simplify_R()
p1.clean_up()
print(p1.to_string())


p2 = ProductTerm()
elements = [R("z"), R("z"),
            Delta("alpha", "beta"),
            MultipoleMoment(["alpha"]),
            MultipoleMoment(["beta"])
            ]
p2.set_elements(elements, prefactor=-1)
p2.simplify_delta()
p2.simplify_R()
p2.clean_up()
print(p2.to_string())

