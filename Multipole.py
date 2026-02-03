import sympy as sp
from sympy import pretty

from Index import Index


class MultipoleMoment:

    def __init__(self, indices:list[str]):

        if len(indices) == 0:
            raise Exception("no Charge")
        if len(indices) == 1:
            self.symbol = sp.Symbol("mu")
        if len(indices) == 2:
            self.symbol = sp.Symbol("Theta")
        if len(indices) == 3:
            self.symbol = sp.Symbol("Omega")
        if len(indices) == 4:
            self.symbol = sp.Symbol("Phi")
        if len(indices) > 4:
            raise Exception("not yet implemented")

        self.indices = [Index(i) for i in indices]


    def to_string(self):
        return pretty(self.symbol) + "_" + "".join([pretty(i.index) for i in self.indices])

    def has_index(self, test_index:Index):
        return test_index.index in [i.index for i in self.indices]

    def replace_index(self, to_be_replaced: Index, replacement: Index):
        if not self.has_index(test_index=to_be_replaced):
            return
        new_indices = []
        for i in self.indices:
            if to_be_replaced.index == i.index:
                new_indices.append(replacement)
            else:
                new_indices.append(i)

        new_indices = sorted([str(i.index) for i in new_indices])
        self.indices = [Index(i) for i in new_indices]


    def is_zero(self):
        if len(self.indices) == 1 and str(self.indices[0].index) in ["x", "z"]:
            return True
        return False
