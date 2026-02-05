import random
import unittest

from mathematics_of_terms.add_product_terms import add_product_terms
from term_components.Delta import Delta
from term_components.Multipole import MultipoleMoment
from term_components.ProductTerm import ProductTerm
from term_components.R import R


class TestMathematics(unittest.TestCase):

    p1 = ProductTerm()
    p2 = ProductTerm()


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
        assert len(result) == 0

        self.p1.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=3)
        self.p2.set_elements(elements=[MultipoleMoment(["alpha"])], prefactor=1)
        result = add_product_terms(self.p1, self.p2)
        assert len(result) == 1
        assert result[0].to_string() == "+ 4 * μ_α"