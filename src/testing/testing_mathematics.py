import unittest
from itertools import permutations

from src.MultipoleInteraction import MultipoleInteraction
from src.mathematics_of_terms.add_product_terms import add_product_terms, addable_product_terms
from src.mathematics_of_terms.sum_up_list import sum_up_list, add_fictitious, addable_fictitious
from src.term_components.Delta import Delta
from src.term_components.Multipole import MultipoleMoment
from src.term_components.ProductTerm import ProductTerm
from src.term_components.R import R


class TestMathematics(unittest.TestCase):

    p1 = ProductTerm()
    p2 = ProductTerm()
    p3 = ProductTerm()
    p4 = ProductTerm()

    dummy_collection = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)


    def test_sort_elements(self):
        self.p1.set_elements(elements=[MultipoleMoment(["gamma"]), MultipoleMoment(["alpha", "beta"])], prefactor=1)
        self.p1.sort_elements()
        assert [i.to_string() for i in self.p1.factors] == ['μ_γ', 'Θ_αβ']

        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "gamma"]), MultipoleMoment(["beta"]) ], prefactor=1)
        self.p1.sort_elements()
        assert [i.to_string() for i in self.p1.factors] == ['μ_β', 'Θ_αγ']



    def test_addable_single_components(self):
        self.p1.set_elements(elements=[], prefactor=1)
        self.p2.set_elements(elements=[], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=2)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=0)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=0)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["beta"])], prefactor=2)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["beta", "alpha"])], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "beta", "gamma"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["beta", "alpha"])], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "beta", "gamma"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["beta", "alpha", "gamma"])], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("beta")], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("alpha")], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[Delta("alpha","beta")], prefactor=1)
        self.p2.set_elements(elements=[Delta("alpha", "beta")], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[Delta("alpha", "beta")], prefactor=1)
        self.p2.set_elements(elements=[Delta("alpha", "gamma")], prefactor=-1)
        result = self.p1.addable(self.p2)
        assert result == False


    def test_addable_combining_same_class(self):
        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), MultipoleMoment(["beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"]), MultipoleMoment(["beta"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), MultipoleMoment(["beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"]), MultipoleMoment(["alpha"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), MultipoleMoment(["alpha", "beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"]), MultipoleMoment(["alpha"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[MultipoleMoment(["gamma"]), MultipoleMoment(["alpha", "beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["gamma"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["gamma"]), MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["delta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["delta"]), MultipoleMoment(["gamma"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["beta", "beta", "beta"]), MultipoleMoment(["beta", "beta"]), MultipoleMoment(["beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["beta", "beta"]), MultipoleMoment(["beta", "beta", "beta"]), MultipoleMoment(["beta"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["beta", "beta", "gamma"]), MultipoleMoment(["beta", "beta"]), MultipoleMoment(["beta"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["beta", "beta"]), MultipoleMoment(["beta", "beta", "beta"]), MultipoleMoment(["beta"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        # todo: R and Delta
        self.p1.set_elements(elements=[R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("alpha")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[R("beta"), R("alpha")], prefactor=1)
        self.p2.set_elements(elements=[R("alpha"), R("beta")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[R("beta"), R("alpha")], prefactor=2)
        self.p2.set_elements(elements=[R("alpha"), R("beta"), R("alpha")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        r1 = R("beta")
        r2 = R("beta")
        r1.exponent = 2
        self.p1.set_elements(elements=[r1], prefactor=2)
        self.p2.set_elements(elements=[r2], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[Delta("beta", "alpha")], prefactor=2)
        self.p2.set_elements(elements=[Delta("alpha","beta")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[Delta("beta", "alpha")], prefactor=2)
        self.p2.set_elements(elements=[Delta("beta","beta")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False


    def test_addable_combining_different_classes(self):
        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("alpha"), MultipoleMoment(["beta"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("beta"), MultipoleMoment(["alpha"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == True

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("beta"), MultipoleMoment(["alpha"]), Delta("alpha", "beta")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("beta"), MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False

        r = R("beta")
        r.exponent = 2
        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        self.p2.set_elements(elements=[r, MultipoleMoment(["alpha"])], prefactor=1)
        result = self.p1.addable(self.p2)
        assert result == False



    def test_add_product_terms(self):
        self.p1.set_elements(elements=[], prefactor=1)
        self.p2.set_elements(elements=[], prefactor=1)
        result = add_product_terms(self.p1, self.p2)
        assert len(result) == 1 # non-vanishing so that terms 0 + 0 are addable

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=3)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=1)
        result = add_product_terms(self.p1, self.p2)
        assert len(result) == 1
        assert result[0].to_string() == "+ 4 * μ_α"

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"]), R("beta")], prefactor=1)
        self.p2.set_elements(elements=[R("beta"), MultipoleMoment(["alpha"])], prefactor=1)
        result = add_product_terms(self.p1, self.p2)
        assert len(result) == 1
        assert result[0].to_string() == "+ 2 * R_β * μ_α"

        self.p1.set_elements(elements=[MultipoleMoment(["alpha","beta"]), R("beta"), Delta("beta", "gamma")], prefactor=-1)
        self.p2.set_elements(elements=[R("beta"), Delta("beta", "gamma"), MultipoleMoment(["alpha","beta"])], prefactor=-3)
        result = add_product_terms(self.p1, self.p2)
        assert len(result) == 1
        assert result[0].to_string() == "- 4 * R_β * Θ_αβ * delta(β, γ)"

    def test_sum_up_list_basics(self):
        assert [10] == sum_up_list([1, 2, 3, 4],
                                   adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [10] == sum_up_list([1, 2, 3, 0, 4],
                                   adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [6] == sum_up_list([1, 2, 3],
                                  adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [6, -1] == sum_up_list([1, 2, 3, -1],
                                      adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [3, -1] == sum_up_list([1, 2, -1],
                                      adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [-1, 5] == sum_up_list([-1, 2, 3],
                                      adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [6, -1] == sum_up_list([1, 2, -1, 3],
                                      adding_function=add_fictitious, addable_function=addable_fictitious)
        assert [-1, 5, -1] == sum_up_list([-1, 2, -1, 3],
                                          adding_function=add_fictitious, addable_function=addable_fictitious)

    def test_sum_up_list_product_terms(self):
        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "beta"]), R("beta"), Delta("beta", "gamma")],
                             prefactor=-1)
        self.p2.set_elements(elements=[R("beta"), Delta("beta", "gamma"), MultipoleMoment(["alpha", "beta"])],
                             prefactor=-3)
        self.p3.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=3)
        self.p4.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=1)


        for order in list(permutations([self.p1, self.p2])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 1
            assert result[0].to_string() == "- 4 * R_β * Θ_αβ * delta(β, γ)"

        for order in list(permutations([self.p3, self.p4])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 1
            assert result[0].to_string() == "+ 4 * μ_α"

        for order in list(permutations([self.p1, self.p2, self.p3, self.p4])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 2
            strings = [i.to_string() for i in result]
            assert "- 4 * R_β * Θ_αβ * delta(β, γ)" in strings
            assert "+ 4 * μ_α" in strings

        for order in list(permutations([self.p1, self.p2, self.p3])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 2
            strings = [i.to_string() for i in result]
            assert "- 4 * R_β * Θ_αβ * delta(β, γ)" in strings
            assert "+ 3 * μ_α" in strings


        for order in list(permutations([self.p1, self.p3, self.p4])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 2
            strings = [i.to_string() for i in result]
            assert "- 1 * R_β * Θ_αβ * delta(β, γ)" in strings
            assert "+ 4 * μ_α" in strings


        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "beta"])], prefactor=-1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha", "beta"])], prefactor=-2)
        self.p3.set_elements(elements=[MultipoleMoment(["alpha", "beta"])], prefactor=-3)
        self.p4.set_elements(elements=[MultipoleMoment(["alpha", "beta"])], prefactor=-4)

        for order in list(permutations([self.p2])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 1
            strings = [i.to_string() for i in result]
            assert "- 2 * Θ_αβ" in strings

        for order in list(permutations([self.p1, self.p2, self.p3])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 1
            strings = [i.to_string() for i in result]
            assert "- 6 * Θ_αβ" in strings

        for order in list(permutations([self.p1, self.p2, self.p3, self.p4])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 1
            strings = [i.to_string() for i in result]
            assert "- 10 * Θ_αβ" in strings

        self.p1.set_elements(elements=[MultipoleMoment(["alpha", "beta", "gamma"])], prefactor=1)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha", "beta", "beta"])], prefactor=2)
        self.p3.set_elements(elements=[MultipoleMoment(["alpha", "beta", "beta"])], prefactor=3)
        self.p4.set_elements(elements=[MultipoleMoment(["alpha", "beta", "gamma"])], prefactor=4)

        for order in list(permutations([self.p1, self.p2, self.p3, self.p4])):
            result = sum_up_list(order,
                                 addable_function=addable_product_terms, adding_function=add_product_terms)
            assert len(result) == 2
            strings = [i.to_string() for i in result]
            assert "+ 5 * Ω_αβγ" in strings
            assert "+ 5 * Ω_αββ" in strings

    # def test_add_up(self):
    #     s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
    #     s.simplify()# includes add_up function
    #     assert ("T_αβγδ μ_α Ω_βγδ = R_z^(-7) * [ "
    #             "- 45 * R_z^(4) * μ_y * Ω_yzz "
    #             "+ 3 * R_z^(4) * μ_y * Ω_yδδ "
    #             "+ 6 * R_z^(4) * μ_y * Ω_yββ "
    #             "]") == s.full_string_tensor()



    def test_multiply_out_greek_indices(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "gamma"])], prefactor=1)

        multiply_out_greek_indices(p)


