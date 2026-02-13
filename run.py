from MultipoleInteraction import MultipoleInteraction
from Tensor import Tensor
from latex_helping_functions import save_lines_as_latex_file

##### FINDING MULTIPOLE INTERACTIONS ##################################################################################

for m1 in range(1,5):
    for m2 in range(m1,5):
        try:
            s = MultipoleInteraction(multipole_order_1 = m1, multipole_order_2 = m2)
            s.simplify_in_latex_steps(importable_tex=False)
        except Exception as e:
            print(f"UNBALE TO DERIVE FOR MULTIPOLE RANKS {m1} and {m2} due to {e}", flush=True)
        print(flush=True)


##### TENSOR DEFINITIONS ##############################################################################################
latex_doc = []
for i in range(7):
    t = Tensor(order=i)
    print(t.full_string_tensor_basics())
    latex_doc.append( t.to_tex(importable_tex=True) )

save_lines_as_latex_file(filename="tensor_definitions.tex", importable_tex=False, content=latex_doc)
