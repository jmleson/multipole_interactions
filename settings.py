import sympy as sp
from sympy import pretty


greek_names = [
    "alpha", "beta", "gamma", "delta", "epsilon",
    "zeta", "eta", "theta", "iota", "kappa",
    "lambda", "mu", "nu", "xi", "omicron",
    "pi", "rho", "sigma", "tau", "upsilon",
    # dummy symbols:
    "phi", "chi", "psi", "omega"
]

global_symbols = sp.symbols(" ".join(greek_names))
symbol_order = {}
for key, value in {name: i for i, name in enumerate(greek_names)}.items():
    symbol_order[key] = value
for key, value in {pretty(name): i for i, name in enumerate(global_symbols)}.items():
    symbol_order[key] = value
for key, value in {"x": -3, "y": -2, "z": -1}.items():
    symbol_order[key] = value

print(symbol_order)
