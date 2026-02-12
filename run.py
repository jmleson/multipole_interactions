from SummandTerm import MultipoleInteraction

# s = MultipoleInteraction(multipole_order_1 = 3, multipole_order_2 = 3)
# s.simplify_in_latex_steps()

for m1 in range(1,5):
    for m2 in range(m1,5):
        try:
            s = MultipoleInteraction(multipole_order_1 = m1, multipole_order_2 = m2)
            s.simplify_in_latex_steps(importable_tex=False)
        except Exception as e:
            print(f"UNBALE TO DERIVE FOR MULTIPOLE RANKS {m1} and {m2} due to {e}", flush=True)
        print(flush=True)
