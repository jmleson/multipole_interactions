from Delta import Delta
from ProductTerm import ProductTerm
from R import R
from settings import global_symbols

import itertools

def variations_between_R_and_delta(amount_of_r: int, amount_of_delta: int):
    order = amount_of_r + 2 * amount_of_delta
    symbols = global_symbols[:order]

    results = []

    # 1️⃣ wähle R-Indizes
    for R_indices in itertools.combinations(symbols, amount_of_r):
        rs = [R(str(i)) for i in R_indices]
        remaining = [s for s in symbols if s not in R_indices]

        # 2️⃣ wähle Delta-Indizes
        for delta_flat in itertools.combinations(remaining, 2 * amount_of_delta):

            deltas = []
            # 3️⃣ in Paare aufteilen (sortiert)
            for i in range(amount_of_delta):
                ds = delta_flat[2*i : 2*i + 2]
                deltas.append(
                    Delta(str(ds[0]), str(ds[1]))
                )

            p = ProductTerm()
            p.set_elements(elements=rs + deltas, prefactor=1)
            results.append(p)

    return results



if __name__ == "__main__":
    res = variations_between_R_and_delta(amount_of_r=1, amount_of_delta=1)
    assert len(res) == 3

