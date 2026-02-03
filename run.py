from Delta import Delta
from Multipole import MultipoleMoment
from ProductTerm import ProductTerm
from R import R


m1 = MultipoleMoment(["alpha"])
m2 = MultipoleMoment(["beta"])

p1 = ProductTerm()
elements = [
    R("alpha"),
    R("beta"),
    m1,
    m2
]
p1.set_elements(elements, prefactor=3)
print(p1.to_string())

p2 = ProductTerm()
elements = [R("z"), R("z"),
            Delta("alpha", "beta"),
            m1, m2]
p2.set_elements(elements, prefactor=-1)
print(p2.to_string())

