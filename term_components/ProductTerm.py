from term_components.Delta import Delta
from term_components.Multipole import MultipoleMoment
from term_components.R import R


class ProductTerm:

    def __init__(self):
        self.factors = []
        self.prefactor = 1

    def set_elements(self, elements:list, prefactor:float):
        for element in elements:
            if isinstance(element, R):
                self.factors.append(element)
            if isinstance(element, MultipoleMoment):
                self.factors.append(element)
            if isinstance(element, Delta):
                self.factors.append(element)
        self.prefactor = prefactor


    def to_string(self):
        if self.prefactor >= 0:
            strings = [f"+ {self.prefactor}"]
        else:
            strings = [f"- {abs(self.prefactor)}"]
        for factor in self.factors:
            strings.append( factor.to_string() )
        return " * ".join(strings)

    def to_latex(self):
        if self.prefactor >= 0:
            strings = [f"+ {self.prefactor}"]
        else:
            strings = [f"- {abs(self.prefactor)}"]
        for factor in self.factors:
            strings.append( factor.to_latex() )
        return " \cdot ".join(strings)

    def simplify(self, type):
        for x in self.factors:
            if isinstance(x, type):
                to_be_replaced, replacement = x.simplify()

                if to_be_replaced is not None and replacement is not None:
                    for i in self.factors:
                        if i.has_index(to_be_replaced):
                            i.replace_index(to_be_replaced=to_be_replaced, replacement=replacement)

                    to_be_replaced, replacement = x.simplify()

    def simplify_delta(self):
        return self.simplify(type=Delta)
        # for d in self.factors:
        #     if isinstance(d, Delta):
        #         to_be_replaced, replacement = d.simplify()
        #
        #         for i in self.factors:
        #             if i.has_index(to_be_replaced):
        #                 i.replace_index(to_be_replaced=to_be_replaced, replacement=replacement)

    def simplify_R(self):
        return self.simplify(type=R)
        # for r in self.factors:
        #     if isinstance(r, R):
        #         to_be_replaced, replacement = r.simplify()
        #
        #         for i in self.factors:
        #             if i.has_index(to_be_replaced):
        #                 i.replace_index(to_be_replaced=to_be_replaced, replacement=replacement)

    def simplify_Multipole(self):
        return self.simplify(type=MultipoleMoment)

    def clean_up(self, set_to_zero:bool = True):
        new_factors = []
        distances = {}
        for factor in self.factors:
            if set_to_zero and factor.is_zero():
                # print("factor is 0:", factor.to_string(), flush=True)
                self.prefactor = 0
                self.factors = []
                return
            if isinstance(factor, Delta):
                result = factor.evaluate()
                if result is None:
                    new_factors.append(factor)
                else:
                    self.prefactor *= result
            elif isinstance(factor, R):
                if str(factor.index.index) not in distances:
                    distances[str(factor.index.index)] = 0
                distances[str(factor.index.index)] += factor.exponent
            else:
                new_factors.append(factor)

        for index, amount in distances.items():
            r = R(str(index))
            r.exponent = amount
            new_factors.append(r)

        self.factors = new_factors
        return




# def simplify_delta(p:ProductTerm, delta:Delta):
#     identical_values = sorted([delta.index1, delta.index2])
#     to_be_replaced, replacement = identical_values[0], identical_values[1]
#
#     for i in p.factors:
#         if i.has_index(to_be_replaced):
#             i.replace_index(to_be_replaced=to_be_replaced, replacement=replacement)
