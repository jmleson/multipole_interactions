from Index import Index
import sympy as sp
from sympy import pretty

from settings import symbol_order


class Delta:

    def __init__(self, index_sign_1:str, index_sign_2:str):
        indices = sorted([index_sign_1, index_sign_2], key=lambda s: symbol_order[str(s)])
        self.index1 = Index(indices[0])
        self.index2 = Index(indices[1])

    def has_index(self, test_index:Index):
        chance_1 = test_index.index == self.index1.index
        chance_2 = test_index.index == self.index2.index
        return chance_1 or chance_2

    def to_string(self):
        delta = sp.Symbol("delta")
        return str(delta) + "(" + pretty(self.index1.index) + ", " + pretty(self.index2.index) + ")"


    def replace_index(self, to_be_replaced:Index, replacement:Index):
        if not self.has_index(test_index=to_be_replaced):
            return

        chance_1 = to_be_replaced.index == self.index1.index
        if chance_1:
            self.index1 = Index(str(replacement.index))
        chance_2 = to_be_replaced.index == self.index2.index
        if chance_2:
            self.index2 = Index(str(replacement.index))

    def evaluate(self):
        if self.index1.index == self.index2.index:
            return 1
        if self.index1.index in ["x", "y", "z"] and self.index2.index in ["x", "y", "z"]:
            if self.index1.index != self.index2.index:
                return 0
        return None # unknown yet


    def simplify(self):
        identical_values = sorted([self.index1, self.index2])
        to_be_replaced, replacement = identical_values[1], identical_values[0]
        return to_be_replaced, replacement

    def is_zero(self):
        if self.evaluate() == 0:
            return True
        return False



