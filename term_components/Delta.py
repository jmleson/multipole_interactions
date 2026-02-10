
import sympy as sp


from settings import symbol_order
from term_components.Index import Index


class Delta:

    def __init__(self, index_sign_1:str, index_sign_2:str):
        indices = sorted([index_sign_1, index_sign_2], key=lambda s: symbol_order[str(s)])
        self.index1 = Index(indices[0])
        self.index2 = Index(indices[1])

    def get_indices(self):
        return [self.index1, self.index2]

    def has_index(self, test_index:Index):
        chance_1 = test_index == self.index1
        chance_2 = test_index == self.index2
        return chance_1 or chance_2

    def to_string(self):
        delta = sp.Symbol("delta")
        return str(delta) + "(" + self.index1.to_string() + ", " + self.index2.to_string() + ")"

    def to_latex(self):
        delta = r"\delta{}"
        return delta + r"\left(" + self.index1.to_latex() + ", " + self.index2.to_latex() + r"\right)"


    def replace_index(self, to_be_replaced:Index, replacement:Index):
        if not self.has_index(test_index=to_be_replaced):
            return

        chance_1 = to_be_replaced.index == self.index1.index
        if chance_1:
            self.index1 = Index(str(replacement.index))
        chance_2 = to_be_replaced.index == self.index2.index
        if chance_2:
            self.index2 = Index(str(replacement.index))

        # check if sorting still in order:
        indices = sorted([self.index1, self.index2])
        if indices != [self.index1, self.index2]:
            #exchange indices
            self.index1, self.index2 = self.index2, self.index1

    def evaluate(self):
        if self.index1.index == self.index2.index:
            return 1
        if self.index1.index in ["x", "y", "z"] and self.index2.index in ["x", "y", "z"]:
            if self.index1.index != self.index2.index:
                return 0
        return None # unknown yet


    def simplify(self) -> list[dict]:
        coordinates = [i for i in [self.index1, self.index2] if str(i.index) in ["x", "y", "z"]]
        greek_indices = [i for i in [self.index1, self.index2] if str(i.index) not in ["x", "y", "z"]]
        if len(coordinates) >= 2:
            return []
        if len(coordinates) == 1:
            return [{"to_be_replaced": greek_indices[0], "replacement": coordinates[0], "dummy": False}]

        dummy_symbol = Index(str("omega"))
        return [
            {"to_be_replaced": greek_indices[1], "replacement": dummy_symbol, "dummy": True},
            {"to_be_replaced": greek_indices[0], "replacement": dummy_symbol, "dummy": True}
        ]


    def is_zero(self):
        if self.evaluate() == 0:
            return True
        return False

    def __eq__(self, other):
        if not isinstance(other, Delta):
            return NotImplemented
        return self.index1 == other.index1 and self.index2 == other.index2

    def __lt__(self, other):
        # sorting
        if not isinstance(other, Delta):
            return NotImplemented
        return (self.index1, self.index2) < (other.index1, other.index2)



