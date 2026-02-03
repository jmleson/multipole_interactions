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
                string_list = sorted([str(i.index) for i in [index1,index2,index3]])
                index_string = "".join(string_list)
                if index_string not in ["xxy", "yyy", "yzz"]:
                    return True
        elif len(self.indices) == 4:
            # hexadecapole:
            index1 = self.indices[0]
            index2 = self.indices[1]
            index3 = self.indices[2]
            index4 = self.indices[3]
            if index1.is_coordinate() and index2.is_coordinate() and index3.is_coordinate() and index4.is_coordinate():
                string_list = sorted([str(i.index) for i in [index1, index2, index3]])
                index_string = "".join(string_list)
                if index_string not in ["xxxx", "xxyy", "xxzz", "yyyy", "yyzz", "zzzz"]:
                    return True
            pass
        else:
            raise Exception("not yet implemented")
        return False


    def simplify(self):
        to_be_replaced, replacement = None, None
        if len(self.indices) == 1:
            #Dipole
            to_be_replaced, replacement = self.indices[0], Index("y")
            #
        elif len(self.indices) == 2:
            #Quadrupole
            index1 = self.indices[0]
            index2 = self.indices[1]
            if index1.is_coordinate() and not index2.is_coordinate():
                to_be_replaced, replacement = index2, index1
            elif not index1.is_coordinate() and index2.is_coordinate():
                to_be_replaced, replacement = index1, index2
            else:
                existing = sorted([str(index1.index), str(index2.index)])
                to_be_replaced, replacement = existing[1], existing[0]
            #
        elif len(self.indices) == 3:
            pass
            #
        elif len(self.indices) == 4:
            pass
            #
        else:
            raise Exception("not yet implemented")

        return to_be_replaced, replacement