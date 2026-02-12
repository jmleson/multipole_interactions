from SummandTerm import MultipoleInteraction
from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms


for m1 in range(1,5):
    for m2 in range(m1,5):
        if sorted([m1,m2]) in ( [1,3], [2,2] ):
            continue
        print("ranks", m1, m2, flush=True )
        s = MultipoleInteraction(multipole_order_1 = m1, multipole_order_2 = m2)
        s.simplify_in_latex_steps()

        t = MultipoleInteraction(multipole_order_1 = m1, multipole_order_2 = m2)
        t.simplify_long_way()

        assert s == t


# s = MultipoleInteraction(multipole_order_1 = 1, multipole_order_2 = 3)
# s.simplify_in_latex_steps()
# t = MultipoleInteraction(multipole_order_1 = 1, multipole_order_2 = 3)
# t.simplify_long_way()
#
# print(s.full_latex_tensor())
# print(t.full_latex_tensor())