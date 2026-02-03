import copy
import unittest
from Delta import Delta
from Index import Index
from Multipole import MultipoleMoment
from ProductTerm import ProductTerm
from R import R



class TestingBasics(unittest.TestCase):
    m1 = MultipoleMoment(["alpha"])
    m2 = MultipoleMoment(["beta"])

    p1 = ProductTerm()
    elements = [
        R("alpha"), R("beta"),
        m1, m2
    ]
    p1.set_elements(elements, prefactor=3)

    p2 = ProductTerm()
    elements = [R("z"), R("z"),
                Delta("alpha", "beta"),
                m1, m2]
    p2.set_elements(elements, prefactor=-1)


    def test_ProductTerm(self):
        p1, p2 = copy.deepcopy(self.p1), copy.deepcopy(self.p2)
        assert "3 * R_α * R_β * μ_α * μ_β" == p1.to_string()
        assert "-1 * R_z * R_z * delta(α, β) * μ_α * μ_β" == p2.to_string()

        p2.clean_up(set_to_zero=True)
        assert "-1 * delta(α, β) * μ_α * μ_β * R_z^(2)" == p2.to_string()

        p1.factors.append(R("y"))
        p1.clean_up(set_to_zero=True)
        assert "0" == p1.to_string()


    def test_replace_index_R(self):
        r = R("alpha")
        r.replace_index(Index("gamma"), Index("delta"))
        assert r.to_string() == "R_α"

        r = R("alpha")
        r.replace_index(Index("alpha"), Index("gamma"))
        assert r.to_string() == "R_γ"


    def test_replace_index_Delta(self):
        d = Delta("alpha", "beta")
        d.replace_index(Index("gamma"), Index("delta"))
        assert d.to_string() == "delta(α, β)"

        d = Delta("alpha", "beta")
        d.replace_index(Index("alpha"), Index("gamma"))
        assert d.to_string() == "delta(γ, β)"

        d = Delta("alpha", "beta")
        d.replace_index(Index("beta"), Index("gamma"))
        assert d.to_string() == "delta(α, γ)"

    def test_replace_index_Multipole(self):
        m = MultipoleMoment(["alpha"])
        m.replace_index(Index("gamma"), Index("delta"))
        assert m.to_string() == "μ_α"

        m = MultipoleMoment(["alpha"])
        m.replace_index(Index("alpha"), Index("beta"))
        assert m.to_string() == "μ_β"

        m = MultipoleMoment(["alpha", "beta"])
        m.replace_index(Index("alpha"), Index("beta"))
        assert m.to_string() == "Θ_ββ"

        m = MultipoleMoment(["alpha", "beta"])
        m.replace_index(Index("alpha"), Index("gamma"))
        assert m.to_string() == "Θ_βγ"


    def test_SimplifyDelta(self):
        p1, p2 = copy.deepcopy(self.p1), copy.deepcopy(self.p2)

        p1.simplify_delta()
        assert "3 * R_α * R_β * μ_α * μ_β" == p1.to_string()

        p2.simplify_delta()
        assert "-1 * R_z * R_z * delta(β, β) * μ_β * μ_β" == p2.to_string()
        p2.clean_up(set_to_zero=True)
        assert "-1 * μ_β * μ_β * R_z^(2)" == p2.to_string()



    def test_SimplifyR(self):
        p1, p2 = copy.deepcopy(self.p1), copy.deepcopy(self.p2)

        # p1.simplify_R()
        # assert "3 * R_z * R_z * μ_z * μ_z" == p1.to_string()
        # p1.clean_up(set_to_zero=False)
        # assert "3 * μ_z * μ_z * R_z^(2)" == p1.to_string()

        p2.simplify_R()
        assert "-1 * R_z * R_z * delta(α, β) * μ_α * μ_β" == p2.to_string()
        p2.clean_up(set_to_zero=True)
        assert "-1 * delta(α, β) * μ_α * μ_β * R_z^(2)" == p2.to_string()

        p3 = ProductTerm()
        elements = [R("alpha"), R("beta"),
                    Delta("alpha", "beta"),
                    self.m1, self.m2]
        p3.set_elements(elements, prefactor=-2)
        p4 = copy.deepcopy(p3)
        p3.simplify_R()
        assert "-2 * R_z * R_z * delta(z, z) * μ_z * μ_z" == p3.to_string()

        p4.simplify_R()
        p4.clean_up(set_to_zero=False)
        assert "-2 * μ_z * μ_z * R_z^(2)" == p4.to_string()
        p4.clean_up(set_to_zero=True)
        assert p4.to_string() == "0"
