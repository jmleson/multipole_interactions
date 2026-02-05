from Delta import Delta
from ProductTerm import ProductTerm
from R import R
from settings import global_symbols

import itertools


def pairings(lst):
    # alle Paarungen der liste lst
    if not lst:
        yield []
        return

    first = lst[0]
    for i in range(1, len(lst)):
        pair = (first, lst[i])
        for rest in pairings(lst[1:i] + lst[i+1:]):
            yield [pair] + rest


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
            for pairing in pairings(list(delta_flat)):
                deltas = [
                    Delta(str(a), str(b))
                    for a, b in pairing
                ]

                p = ProductTerm()
                p.set_elements(elements=rs + deltas, prefactor=1)
                results.append(p)

    return results





if __name__ == "__main__":
    res = variations_between_R_and_delta(amount_of_r=1, amount_of_delta=1)
    assert len(res) == 3

