
import sympy as sp



class Index:

    def __init__(self, index_sign:str):
        self.index = sp.Symbol(index_sign)

    def __lt__(self, other):
        if not isinstance(other, Index):
            return NotImplemented
        return str(self.index) < str(other.index)

    def is_coordinate(self):
        if str(self.index) in ["x", "y", "z"]:
            return True