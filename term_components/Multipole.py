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

        self.molecule = None


    def get_indices(self):
        return self.indices

    def get_greek_indices(self):
        greek_indices = [str(i.index) for i in self.get_indices()]
        greek_indices = [i for i in greek_indices if i not in {"x", "y", "z"}]
        greek_indices = set(greek_indices)
        greek_indices = list([Index(i) for i in greek_indices])
        return greek_indices

    def to_string(self):
        if self.molecule:
            return pretty(self.symbol) + r"^" + self.molecule + "_" + "".join([i.to_string() for i in self.indices])
        return pretty(self.symbol) + "_" + "".join([i.to_string() for i in self.indices])

    def to_latex(self):
        if self.molecule:
             return r"\ "[0]+str(self.symbol) + "{}"+ r"^{" + self.molecule + r"}"+ r"_{" + "".join([i.to_latex() for i in self.indices]) + "}"
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
                    index2 = Index(indices[1])
                    if not index1.is_coordinate() and index2.is_coordinate():
                        return [index1]
                    if index1.is_coordinate() and not index2.is_coordinate():
                        return [index2]
                    if not index1.is_coordinate() and not index2.is_coordinate():
                        index1.separate_traceless = True
                        index1.separate_traceless = True
                        return [index1, index2]
        return []




    def simplify(self) -> list[dict]:
        replacements = []

        if len(self.indices) == 1:
            if self.indices[0].is_coordinate():
                return replacements
            #Dipole
            replacements.append( {"to_be_replaced": self.indices[0], "replacement": Index("y"), "dummy": False} )

        elif len(self.indices) == 2:
            #Quadrupole
            index1 = self.indices[0]
            index2 = self.indices[1]
            if index1.is_coordinate() and index2.is_coordinate():
                return replacements
            elif index1.is_coordinate() and not index2.is_coordinate():
                replacements.append({"to_be_replaced": index2, "replacement": index1, "dummy": False})
            elif not index1.is_coordinate() and index2.is_coordinate():
                replacements.append({"to_be_replaced": index1, "replacement": index2, "dummy": False})
            else:
                # existing = sorted([str(index1.index), str(index2.index)], key=lambda s: symbol_order[str(s)])
                dummy = Index("omega")
                replacements.append({"to_be_replaced": Index(str(index1.index)), "replacement": dummy, "dummy": True})
                replacements.append({"to_be_replaced": Index(str(index2.index)), "replacement": dummy, "dummy": True})
            #
        elif len(self.indices) == 3:
            coordinates = [i for i in self.indices if i.is_coordinate()]
            greek_symbols = sorted([i for i in self.indices if not i.is_coordinate()])
            if len(coordinates) == 3:
                return replacements
            if len(coordinates) == 2 and len(greek_symbols) == 1:
                if Index("y") in coordinates:
                    coordinates.remove(Index("y"))# removes 1st y from list
                    replacements.append({"to_be_replaced": greek_symbols[0],
                                         "replacement": coordinates[0], "dummy": False })
                else:
                    replacements.append({"to_be_replaced": greek_symbols[0], "replacement": Index("y"), "dummy": False})
            if len(coordinates) == 1 and len(greek_symbols) == 2:
                if coordinates[0] == Index("y") and len(set([str(i.index) for i in greek_symbols])) == 2:
                    dummy = Index("omega")
                    replacements.append({"to_be_replaced": greek_symbols[1], "replacement": dummy, "dummy": True})
                    replacements.append({"to_be_replaced": greek_symbols[0], "replacement": dummy, "dummy": True})
            if len(coordinates) == 0 and len(greek_symbols) == 3:
                counts = Counter(i.index for i in greek_symbols)
                if sorted(counts.values()) == [1, 2]:
                    single_index = next((k for k, v in counts.items() if v == 1), None)
                    if single_index is not None:
                        replacements.append({"to_be_replaced": Index(str(single_index)), "replacement": Index("y"), "dummy": False})
            #
        elif len(self.indices) == 4:
            coordinates = [i for i in self.indices if i.is_coordinate()]
            greek_symbols = sorted([i for i in self.indices if not i.is_coordinate()])

            if len(coordinates) == 4:
                return replacements
            if len(coordinates) == 3:
                counts = Counter(i.index for i in coordinates)
                if sorted(counts.values()) == [1, 2]:
                    single_index = next((k for k, v in counts.items() if v == 1), None)
                    if single_index is not None:
                        replacements.append(
                            {"to_be_replaced": greek_symbols[0], "replacement": Index(str(single_index)), "dummy": False})

            if len(coordinates) == 2 and len(greek_symbols) == 2:
                if coordinates[0] != coordinates[1]:
                    pass
                    # # greek symbols have to be identical to the coordinates
                    # replacements.append({"to_be_replaced": greek_symbols[0], "replacement": coordinates[0]})
                    # replacements.append({"to_be_replaced": greek_symbols[1], "replacement": coordinates[1]})
                    # # OR:
                    # replacements.append({"to_be_replaced": greek_symbols[0], "replacement": coordinates[1]})
                    # replacements.append({"to_be_replaced": greek_symbols[1], "replacement": coordinates[0]})
                    # # -> leads to factor 2 (2 summand terms remain)
                else:
                    if len(set([str(i.index) for i in greek_symbols])) == 2:
                        # otherwise the indices are already indicating the same value
                        dummy = Index("omega")
                        replacements.append({"to_be_replaced": greek_symbols[1], "replacement": dummy, "dummy": True})
                        replacements.append({"to_be_replaced": greek_symbols[0], "replacement": dummy, "dummy": True})
            if len(coordinates) == 1:
                counts = Counter(i.index for i in greek_symbols)
                if sorted(counts.values()) == [1, 2]:
                    single_index = next((k for k, v in counts.items() if v == 1), None)
                    if single_index is not None:
                        replacements.append(
                            {"to_be_replaced": Index(str(single_index)), "replacement": coordinates[0], "dummy": False})
            if len(coordinates) == 0:
                counts = Counter(i.index for i in greek_symbols)
                if sorted(counts.values()) == [1, 1, 2]:
                    double_index = next((k for k, v in counts.items() if v == 2), None)
                    if double_index is not None:
                        dummy = Index("omega")
                        single_indices = [i for i in greek_symbols if str(i.index) != str(double_index)]
                        for single_index in single_indices:
                            replacements.append({"to_be_replaced": single_index, "replacement": dummy, "dummy": True})
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
        if self.molecule and not other.molecule or not self.molecule and other.molecule:
            raise Exception("not possible")
        return self.indices == other.indices and self.prefactor_in_potential == other.prefactor_in_potential and self.molecule == other.molecule


    def __lt__(self, other):
        # sorting
        if not isinstance(other, MultipoleMoment):
            return NotImplemented
        if self.molecule and not other.molecule or not self.molecule and other.molecule:
            raise Exception("unable to compare")
        return (self.molecule, len(self.indices), self.indices) < (other.molecule, len(other.indices), other.indices)