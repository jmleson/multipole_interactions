import sympy as sp


from settings import greek_names, symbol_order
from term_components.Delta import Delta
from term_components.Index import Index
from term_components.Multipole import MultipoleMoment
from term_components.R import R


class ProductTerm:

    def __init__(self):
        self.factors = []
        self.prefactor = 1

    def set_elements(self, elements:list, prefactor:float):
        self.factors = []
        for element in elements:
            if isinstance(element, R):
                self.factors.append(element)
            if isinstance(element, MultipoleMoment):
                self.factors.append(element)
            if isinstance(element, Delta):
                self.factors.append(element)
        self.prefactor = prefactor
        if len(elements) == 0:
            self.prefactor = 0

    def get_indices(self):
        indices = []
        for i in self.factors:
            indices.extend(i.get_indices())
        return indices


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
            strings = [f"+ {sp.latex(self.prefactor)}"]
        else:
            strings = [f"- {sp.latex(abs(self.prefactor))}"]
        for factor in self.factors:
            strings.append( factor.to_latex() )
        return r" \cdot ".join(strings)

    def simplify(self, type: R|Delta|MultipoleMoment, max_order:int=None):
        """
        :param type: object class (of element in self.factors) that should be simplified
        :param max_order: for multipoles only -> cut simplification if order of multipole higher than max_order
        :return:
        """
        for x in self.factors:
            if isinstance(x, type):
                if type == MultipoleMoment and max_order is not None and len(x.indices) > max_order:
                    continue

                possible_replacements = x.simplify()
                self._replace_given_replacement(possible_replacements=possible_replacements)


    def _replace_given_replacement(self, possible_replacements:list[dict]):
        # checking dummy set-up:
        to_be_replaced_dummies = []
        for replacement in possible_replacements:
            if replacement["dummy"]:
                if replacement["replacement"] not in to_be_replaced_dummies:
                    to_be_replaced_dummies.append(replacement["replacement"])
                for i in self.factors:
                    if i.has_index(replacement["replacement"]):
                        raise Exception("dummy should be index that does not occur otherwise (that's the whole point!)")

        # actual replacement
        for replacement in possible_replacements:
            if replacement["to_be_replaced"] is not None and replacement["replacement"] is not None:
                for i in self.factors:
                    if i.has_index(replacement["to_be_replaced"]):
                        i.replace_index(to_be_replaced=replacement["to_be_replaced"],
                                        replacement=replacement["replacement"])

        # print("\tafter dummy:", self.to_string(), flush=True)#TODO remove

        # cleaning up dummies:
        if len(to_be_replaced_dummies) > 0:
            indices = [str(i.index) for i in self.get_indices()]
            first_not_used_greek_index = [i for i in greek_names if i not in indices]
            if len(first_not_used_greek_index):
                for d in range(len(to_be_replaced_dummies)):
                    dummy = to_be_replaced_dummies[d]
                    for i in self.factors:
                        i.replace_index(to_be_replaced=dummy, replacement=Index(first_not_used_greek_index[d]))

    def simplify_delta(self):
        return self.simplify(type=Delta)

    def simplify_R(self):
        return self.simplify(type=R)

    def simplify_Multipole(self, max_order:int=None):
        return self.simplify(type=MultipoleMoment, max_order=max_order)

    def reduce_indices(self, consider_index_sorting:bool):
        """
        moves terms of higher indices to lower ones, e.g. mu_beta -> mu_alpha
        (! not the smartest function, it does only re-name indices when their position in the alphabet is higher than the number of used indices )
        :return:
        """
        self.sort_elements()# needed to keep replacements independent of initial order of elements

        used_indices = []
        for i in self.factors:
            used_indices.extend(i.get_indices())
        used_indices = sorted(used_indices)
        included_coordinates = [str(i.index) for i in used_indices if i.is_coordinate()]
        included_coordinates = [Index(i) for i in set(included_coordinates)]
        used_indices = [str(i.index) for i in used_indices]

        used_indices = list(dict.fromkeys(used_indices))
        index_map = {}
        greek_iter = iter(greek_names)
        for old_index in used_indices:
            # Prüfen, ob der alte Index schon in den kleinsten verfügbaren ist
            if old_index in greek_names[:len(used_indices)-len(included_coordinates)]:
                # Index ist schon optimal, nichts tun
                continue
            if Index(old_index).is_coordinate():
                continue
            # Bndernfalls auf den nächsten freien griechischen Namen abbilden
            while True:
                new_index = next(greek_iter)
                if new_index not in used_indices and new_index not in index_map.values():
                    index_map[old_index] = new_index
                    break

        for to_be_replaced, replacement in index_map.items():
            to_be_replaced = Index(to_be_replaced)
            replacement = Index(replacement)
            for i in self.factors:
                if i.has_index(replacement):
                    raise Exception("index already exists, by replacement a change in meaning will be produced ")
                if i.has_index(to_be_replaced):
                    i.replace_index(to_be_replaced=to_be_replaced, replacement=replacement)

        if consider_index_sorting:
            self.reduce_indices_sorting()
        return

    def reduce_indices_sorting(self):
        total = self.to_string()
        index_order_ist = []
        for i in total:
            if i in symbol_order.keys() and i not in index_order_ist and i not in ["x", "y", "z"]:
                index_order_ist.append(i)

        index_order_shall = []
        for i in sorted(self.get_indices()):
            if i.to_string() not in index_order_shall and i.to_string() not in ["x", "y", "z"]:
                index_order_shall.append(i.to_string())

        if index_order_shall != index_order_ist:
            possible_replacements = []
            for i in range(len(index_order_shall)):
                replacement = {"to_be_replaced": Index(index_order_shall[-i-1]), "replacement": Index(greek_names[-1-i]), "dummy": True}
                possible_replacements.append(replacement)

            self._replace_given_replacement(possible_replacements)
        if total != self.to_string():
            print("\t\t! function 'reduce_indices_sorting' changed", total , "into", self.to_string(), flush=True)
        return



    def use_traceless_conditions(self):
        for x in range(len(self.factors)):
            term = self.factors[x]
            if not isinstance(term, MultipoleMoment):
                continue

            traceless = term.traceless_condition()
            if len(traceless) == 0:
                continue

            indices = []
            for y in range(len(self.factors)):
                if y != x: #! important to compare indices, because the same multipole may occur twice/...
                    indices.extend(self.factors[y].get_indices())

            index_does_occur = []
            for possible_zero_component in traceless:
                if possible_zero_component.is_coordinate():
                    continue

                if possible_zero_component.separate_traceless:
                    if not possible_zero_component in indices:
                        index_does_occur.append(True)
                else:
                    if possible_zero_component in indices:
                        index_does_occur.append(False)
            if not False in index_does_occur or True in index_does_occur:
                self.prefactor = 0
                self.factors = []
                return


    def includes_duplicates(self) -> tuple[bool,list[Index]]:
        duplicates_included = False
        greek_indices = []
        for i in range(len(self.factors)):
            if not isinstance(self.factors[i], MultipoleMoment):
                continue
            for j in range(i+1,len(self.factors)):
                if self.factors[i] == self.factors[j]:
                    duplicates_included = True
                    greek_indices = self.factors[i].get_greek_indices()

        return duplicates_included, greek_indices

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
        self.sort_elements()
        return

    def sort_elements(self):
        type_order = {R: 0, MultipoleMoment: 1, Delta: 2}
        self.factors = sorted(
            self.factors,
            key=lambda x: (type_order.get(type(x), 99), x)
        )

    def addable(self, other):
        if not isinstance(other, ProductTerm):
            return False
        if len(self.factors) != len(other.factors):
            return False

        self.sort_elements()
        other.sort_elements()
        for i in range(len(self.factors)):
            if type(self.factors[i]) != type(other.factors[i]):
                return False
            if type(self.factors[i]) == MultipoleMoment:
                if self.factors[i].indices != other.factors[i].indices:
                    # print([i.to_string() for i in self.factors],
                    #       [i.to_string() for i in other.factors], flush=True )
                    return False
            elif self.factors[i] != other.factors[i]:
                return False
        return True

