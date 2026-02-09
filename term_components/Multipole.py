import sympy as sp
from sympy import pretty


from settings import symbol_order
from term_components.Index import Index
from collections import Counter

class MultipoleMoment:

    def __init__(self, indices:list[str]):
        if not isinstance(indices, list):
            raise Exception("wrong format of input")

        if len(indices) == 0:
            raise Exception("no Charge")
        elif len(indices) == 1:
            self.symbol = sp.Symbol("mu")
        elif len(indices) == 2:
            self.symbol = sp.Symbol("Theta")
        elif len(indices) == 3:
            self.symbol = sp.Symbol("Omega")
        elif len(indices) == 4:
            self.symbol = sp.Symbol("Phi")
        elif len(indices) > 4:
            raise Exception("not yet implemented")

        self.indices = sorted([Index(i) for i in indices])

        self.prefactor_in_potential = (-1) ** len(self.indices)


    def get_indices(self):
        return self.indices

    def get_greek_indices(self):
        greek_indices = [str(i.index) for i in self.get_indices()]
        greek_indices = [i for i in greek_indices if i not in {"x", "y", "z"}]
        greek_indices = set(greek_indices)
        greek_indices = list([Index(i) for i in greek_indices])
        return greek_indices

    def to_string(self):
        return pretty(self.symbol) + "_" + "".join([i.to_string() for i in self.indices])

    def to_latex(self):
        return r"\ "[0]+str(self.symbol) + "{}_{" + "".join([i.to_latex() for i in self.indices]) + "}"

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

        new_indices = sorted([str(i.index) for i in new_indices], key=lambda s: symbol_order[str(s)])
        self.indices = [Index(i) for i in new_indices]


    def is_zero(self):
        if len(self.indices) == 1:
            # dipole: index must be y
            if str(self.indices[0].index) in ["x", "z"]:
                return True
        elif len(self.indices) == 2:
            # Quadrupole: indices must be identical
            index1 = self.indices[0]
            index2 = self.indices[1]
            if index1.is_coordinate() and index2.is_coordinate() and str(index1.index) != str(index2.index):
                return True
        elif len(self.indices) == 3:
            # octopole:
            index1 = self.indices[0]
            index2 = self.indices[1]
            index3 = self.indices[2]
            if index1.is_coordinate() and index2.is_coordinate() and index3.is_coordinate():
                string_list = sorted([str(i.index) for i in [index1,index2,index3]], key=lambda s: symbol_order[str(s)])
                index_string = "".join(string_list)
                if index_string not in ["xxy", "yyy", "yzz"]:
                    return True
            coordinates = [i for i in [index1, index2, index3] if i.is_coordinate()]
            greek_symbols = [i for i in [index1, index2, index3] if not i.is_coordinate()]
            if len(coordinates) == 1 and Index("y") not in coordinates:
                # can only become xxy / yzz -> at least one y has to be added
                if greek_symbols[0] == greek_symbols[1]:
                    return True

        elif len(self.indices) == 4:
            # hexadecapole:
            index1 = self.indices[0]
            index2 = self.indices[1]
            index3 = self.indices[2]
            index4 = self.indices[3]
            coordinates = [i for i in [index1, index2, index3, index4] if i.is_coordinate()]
            greek_symbols = [i for i in [index1, index2, index3, index4] if not i.is_coordinate()]
            if len(coordinates) == 4: #index1.is_coordinate() and index2.is_coordinate() and index3.is_coordinate() and index4.is_coordinate():
                string_list = sorted([str(i.index) for i in [index1, index2, index3, index4]], key=lambda s: symbol_order[str(s)])
                index_string = "".join(string_list)
                if index_string not in ["xxxx", "xxyy", "xxzz", "yyyy", "yyzz", "zzzz"]:
                    return True
            if len(coordinates) == 2:
                if (coordinates[0] != coordinates[1]  # other indices have to be identical to these -> other indices cannot be identical to eachother
                            and greek_symbols[0] == greek_symbols[1]):
                    return True
        else:
            raise Exception("not yet implemented")
        return False

    def traceless_condition(self) -> list[Index]:
        """
        checks whether a traceless condition can be applied ( len(result) > 0 )
        :return: list of affected indices (in case a traceless condition is fullfilled)
        """
        if len(self.indices) == 2:
            if self.indices[0] == self.indices[1] and not self.indices[0].is_coordinate() and not self.indices[1].is_coordinate():
                return [Index(str(self.indices[0].index))]
        if len(self.indices) == 3:
            # = 0, if xαα / yαα / zαα / or 2x other greek symbol :
            counts = Counter(i.index for i in self.indices)
            if len(counts) == 2:
                values = list(counts.values())
                if sorted(values) == [1, 2]:
                    double_index = None
                    for k, v in counts.items():
                        if v == 2:
                            double_index = k
                            break
                    if double_index is not None and not Index(str(double_index)).is_coordinate():
                        return [Index(str(double_index))]# single_index can be used somewhere else -> still 0 if double_index double
        if len(self.indices) == 4:
            counts = Counter(i.index for i in self.indices)
            if len(counts) == 4: # 4 different
                return []

            if len(counts) == 2:
                values = list(counts.values())
                if values == [2, 2]:
                    indices = [str(i) for i in counts.keys()]
                    index1 = Index(indices[0])
                    index1.separate_traceless = True
                    index2 = Index(indices[1])
                    index1.separate_traceless = True
                    if not index1.is_coordinate() and index2.is_coordinate():
                        return [index1]
                    if index1.is_coordinate() and not index2.is_coordinate():
                        return [index2]
                    if not index1.is_coordinate() and not index2.is_coordinate():
                        return [index1, index2]

        return []




    def simplify(self) -> list[dict]:
        replacements = []
        to_be_replaced, replacement = None, None
        if len(self.indices) == 1:
            if self.indices[0].is_coordinate():
                return replacements
            #Dipole
            replacements.append( {"to_be_replaced": self.indices[0], "replacement": Index("y")} )

        elif len(self.indices) == 2:
            #Quadrupole
            index1 = self.indices[0]
            index2 = self.indices[1]
            if index1.is_coordinate() and index2.is_coordinate():
                return replacements
            elif index1.is_coordinate() and not index2.is_coordinate():
                replacements.append({"to_be_replaced": index2, "replacement": index1})
            elif not index1.is_coordinate() and index2.is_coordinate():
                replacements.append({"to_be_replaced": index1, "replacement": index2})
            else:
                existing = sorted([str(index1.index), str(index2.index)], key=lambda s: symbol_order[str(s)])
                replacements.append({"to_be_replaced": Index(existing[1]), "replacement": Index(existing[0])})
            #
        elif len(self.indices) == 3:
            index1 = self.indices[0]
            index2 = self.indices[1]
            index3 = self.indices[2]
            if index1.is_coordinate() and index2.is_coordinate() and index3.is_coordinate():
                return replacements
            pass
            #
        elif len(self.indices) == 4:
            index1 = self.indices[0]
            index2 = self.indices[1]
            index3 = self.indices[2]
            index4 = self.indices[3]
            if index1.is_coordinate() and index2.is_coordinate() and index3.is_coordinate() and index4.is_coordinate():
                return replacements
            coordinates = [i for i in [index1, index2, index3, index4] if i.is_coordinate()]
            greek_symbols = sorted([i for i in [index1, index2, index3, index4] if not i.is_coordinate()])
            if len(coordinates) == 2:
                if coordinates[0] != coordinates[1]:
                    # greek symbols have to be identical to the coordinates
                    replacements.append({"to_be_replaced": greek_symbols[0], "replacement": coordinates[0]})
                    replacements.append({"to_be_replaced": greek_symbols[1], "replacement": coordinates[1]})
                else:
                    replacements.append({"to_be_replaced": greek_symbols[1], "replacement": greek_symbols[0]})
        else:
            raise Exception("not yet implemented")

        for replacement in replacements:
            if replacement["to_be_replaced"] is not None and not isinstance(replacement["to_be_replaced"], Index):
                raise TypeError("test_index must be an Index")
            if replacement["replacement"] is not None and not isinstance(replacement["replacement"], Index):
                raise TypeError("test_index must be an Index")
        return replacements


    def __eq__(self, other):
        if not isinstance(other, MultipoleMoment):
            return NotImplemented
        return self.indices == other.indices and self.prefactor_in_potential == other.prefactor_in_potential



    def __lt__(self, other):
        # sorting
        if not isinstance(other, MultipoleMoment):
            return NotImplemented
        return (len(self.indices), self.indices) < (len(other.indices), other.indices)