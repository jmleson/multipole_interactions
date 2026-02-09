from term_components.Index import Index


class R:

    def __init__(self, index_sign:str):
        self.index = Index(index_sign)
        self.exponent = 1

    def get_indices(self):
        return [self.index]

    def to_string(self):
        if self.exponent != 1:
            return  "R_" + self.index.to_string() + f"^({self.exponent})"
        return "R_" + self.index.to_string()

    def to_latex(self):
        if self.exponent != 1:
            return "R_{" + self.index.to_latex() + r"}" + r"^{" + str(self.exponent) + r"}"
        return "R_{" + self.index.to_latex() + r"}"

    def has_index(self, test_index:Index):
        if not isinstance(test_index, Index):
            raise TypeError("test_index must be an Index")
        return test_index.index == self.index.index

    def replace_index(self, to_be_replaced: Index, replacement: Index):
        if not self.has_index(test_index=to_be_replaced):
            return
        self.index = Index(str(replacement.index))

    def simplify(self) -> list[dict]:
        to_be_replaced, replacement = self.index, Index("z")
        return [{"to_be_replaced": to_be_replaced, "replacement": replacement}]

    def is_zero(self):
        if str(self.index.index) in ["x", "y"]:
            return True
        return False


    def __eq__(self, other):
        if not isinstance(other, R):
            return NotImplemented
        return self.index == other.index and self.exponent == other.exponent


    def __lt__(self, other):
        # sorting
        if not isinstance(other, R):
            return NotImplemented
        return (self.index, self.exponent) < (other.index, other.exponent)