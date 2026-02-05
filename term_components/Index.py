
import sympy as sp

from settings import symbol_order
from sympy import pretty

class Index:

    def __init__(self, index_sign:str):
        self.index = sp.Symbol(index_sign)

    def __lt__(self, other):
        if not isinstance(other, Index):
            return NotImplemented
        return (
            symbol_order.get(self.index, float("inf")),
            str(self.index),
        ) < (
            symbol_order.get(other.index, float("inf")),
            str(other.index),
        )

    def is_coordinate(self):
        if str(self.index) in ["x", "y", "z"]:
            return True

    def to_string(self):
        return pretty(self.index)

    def to_latex(self):
        if str(self.index) in ["x", "y", "z"]:
            return str(self.index)
        return r"\ "[0] + str(self.index) + r"{}"