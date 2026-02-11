from build_up_tensor_equations.construct_full_tensor import prefactor_expansion

for i in range(6):
    print(i, "->", prefactor_expansion(i) )