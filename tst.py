from SummandTerm import MultipoleInteraction
from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms

construct_all_tensor_terms(order=6)

s = MultipoleInteraction(multipole_order_1 = 3, multipole_order_2 = 3)
s.simplify_in_latex_steps()