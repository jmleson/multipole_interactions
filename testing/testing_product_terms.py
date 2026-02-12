import copy
import unittest
from itertools import permutations

from term_components.Delta import Delta
from term_components.Index import Index
from term_components.Multipole import MultipoleMoment
from term_components.ProductTerm import ProductTerm
from term_components.R import R



class TestingProductTerms(unittest.TestCase):
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
        assert "+ 3 * R_α * R_β * μ_α * μ_β" == p1.to_string()
        assert "- 1 * R_z * R_z * delta(α, β) * μ_α * μ_β" == p2.to_string()

        p2.clean_up(set_to_zero=True)
        assert "- 1 * R_z^(2) * μ_α * μ_β * delta(α, β)" == p2.to_string()

        p1.factors.append(R("y"))
        p1.clean_up(set_to_zero=True)
        assert "+ 0" == p1.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "x", "x"])], prefactor=3)
        p.clean_up()
        assert "+ 0" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "z", "x"])], prefactor=3)
        p.clean_up()
        assert "+ 0" == p.to_string()


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
        assert d.to_string() == "delta(β, γ)"

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

        m = MultipoleMoment(["alpha", "beta", "gamma", "delta"])
        m.replace_index(Index("alpha"), Index("gamma"))
        assert m.to_string() == "Φ_βγγδ"

    def test_is_zero_Oktopole(self):
        m = MultipoleMoment(["x", "y", "z"])
        assert m.is_zero() is True

        m = MultipoleMoment(["alpha", "beta", "gamma"])
        assert m.is_zero() is False

        m = MultipoleMoment(["alpha", "beta", "beta"])
        assert m.is_zero() is False # only 0 because of traceless condition

        m = MultipoleMoment(["x", "beta", "beta"])
        assert m.is_zero() is True
        m = MultipoleMoment(["z", "beta", "beta"])
        assert m.is_zero() is True
        m = MultipoleMoment(["y", "beta", "beta"])
        assert m.is_zero() is False

        m = MultipoleMoment(["y", "alpha", "beta"])
        assert m.is_zero() is False# alpha has to be equal to beta, but that's possible
        m = MultipoleMoment(["x", "alpha", "beta"])
        assert m.is_zero() is False
        m = MultipoleMoment(["z", "alpha", "beta"])
        assert m.is_zero() is False

        m = MultipoleMoment(["y", "y", "beta"])
        assert m.is_zero() is False# beta = y
        m = MultipoleMoment(["x", "x", "beta"])
        assert m.is_zero() is False# beta = y
        m = MultipoleMoment(["z", "z", "beta"])
        assert m.is_zero() is False# beta = y

        m = MultipoleMoment(["y", "x", "beta"])
        assert m.is_zero() is False  # beta = x
        m = MultipoleMoment(["x", "z", "beta"])
        assert m.is_zero() is True
        m = MultipoleMoment(["y", "z", "beta"])
        assert m.is_zero() is False  # beta = z



    def test_SimplifyDelta(self):
        p1, p2 = copy.deepcopy(self.p1), copy.deepcopy(self.p2)

        # no delta included -> no change:
        p1.sort_elements()
        assert p1.to_string() == "+ 3 * R_α * R_β * μ_α * μ_β"
        p1.simplify_delta()
        assert "+ 3 * R_α * R_β * μ_α * μ_β" == p1.to_string()

        # 1 delta included
        p2.sort_elements()
        assert p2.to_string() == "- 1 * R_z * R_z * μ_α * μ_β * delta(α, β)"
        p2.simplify_delta()
        # 1st step: "- 1 * R_z * R_z * μ_ω * μ_ω * delta(ω, ω)" == p2.to_string()
        assert "- 1 * R_z * R_z * μ_α * μ_α * delta(α, α)" == p2.to_string()# 2nd step
        p2.clean_up(set_to_zero=True)
        assert "- 1 * R_z^(2) * μ_α * μ_α" == p2.to_string()

        p = ProductTerm()
        p.set_elements([Delta("alpha", "beta"), R("gamma")], prefactor=-2)
        p.sort_elements()
        assert p.to_string() == "- 2 * R_γ * delta(α, β)"
        p.simplify_delta()
        assert p.to_string() == "- 2 * R_γ * delta(α, α)"

        p = ProductTerm()
        p.set_elements([Delta("x", "beta"), R("beta")], prefactor=-2)
        p.sort_elements()
        assert p.to_string() == "- 2 * R_β * delta(x, β)"
        p.simplify_delta()
        assert p.to_string() == "- 2 * R_x * delta(x, x)"

        p = ProductTerm()
        p.set_elements([Delta("x", "x"), R("x")], prefactor=-2)
        p.sort_elements()
        assert p.to_string() == "- 2 * R_x * delta(x, x)"
        p.simplify_delta()
        assert p.to_string() == "- 2 * R_x * delta(x, x)"

    def test_SimplifyR(self):
        p1, p2 = copy.deepcopy(self.p1), copy.deepcopy(self.p2)

        p2.simplify_R()
        assert "- 1 * R_z * R_z * delta(α, β) * μ_α * μ_β" == p2.to_string()
        p2.clean_up(set_to_zero=True)
        p2.sort_elements()
        assert "- 1 * R_z^(2) * μ_α * μ_β * delta(α, β)" == p2.to_string()

        p3 = ProductTerm()
        elements = [R("alpha"), R("beta"),
                    Delta("alpha", "beta"),
                    self.m1, self.m2]
        p3.set_elements(elements, prefactor=-2)
        p4 = copy.deepcopy(p3)
        p3.simplify_R()
        assert "- 2 * R_z * R_z * delta(z, z) * μ_z * μ_z" == p3.to_string()

        p4.simplify_R()
        p4.clean_up(set_to_zero=False)
        assert "- 2 * R_z^(2) * μ_z * μ_z" == p4.to_string()
        p4.clean_up(set_to_zero=True)
        assert p4.to_string() == "+ 0"


    def test_SimplifyDipoles(self):
        p1, p2 = copy.deepcopy(self.p1), copy.deepcopy(self.p2)

        p1.simplify_Multipole()
        assert "+ 3 * R_y * R_y * μ_y * μ_y" == p1.to_string()

        p2.simplify_Multipole()
        assert "- 1 * R_z * R_z * delta(y, y) * μ_y * μ_y" == p2.to_string()
        p2.clean_up(set_to_zero=True)
        assert "- 1 * R_z^(2) * μ_y * μ_y" == p2.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * μ_y" == p.to_string()
        p.clean_up(set_to_zero=True)
        assert "+ 2 * μ_y" == p.to_string()

        for symbol in ["x", "z"]:
            p = ProductTerm()
            p.set_elements([MultipoleMoment([symbol])], prefactor=2)
            p.simplify_Multipole()
            assert f"+ 2 * μ_{symbol}" == p.to_string()
            p.clean_up(set_to_zero=True)
            assert "+ 0" == p.to_string()

    def test_SimplifyQuadrupoles(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Θ_αα" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "beta"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Θ_xx" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "y"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Θ_yy" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "y"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Θ_xy" == p.to_string()
        p.clean_up(set_to_zero=True)
        assert "+ 0" == p.to_string()

    def test_SimplifyOctopoles(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "y", "alpha"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Ω_xxy" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y", "y", "beta"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Ω_yyy" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "x", "beta"])], prefactor=2)
        p.simplify_Multipole()
        assert "+ 2 * Ω_xxy" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "y"])], prefactor=2)
        assert "+ 2 * Ω_yαβ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * Ω_yαα" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "z"])], prefactor=2)
        assert "+ 2 * Ω_zαβ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * Ω_zαβ" == p.to_string()# two choices remain: alpha = z and beta = y OR alpha = y and beta = z -> factor 2

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "beta"])], prefactor=2)
        assert "+ 2 * Ω_αββ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * Ω_yββ" == p.to_string()

    def test_SimplifyHexadecapoles(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y", "x", "delta", "x"])], prefactor=3)
        p.sort_elements()
        assert "+ 3 * Φ_xxyδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 3 * Φ_xxyy" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "x", "delta", "delta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * Φ_xxδδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * Φ_xxδδ" == p.to_string()  # cannot be replaced further
        p.clean_up(set_to_zero=True)  # assert not p.factors[0].is_zero()
        assert "+ 2 * Φ_xxδδ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "x", "delta", "beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * Φ_xxβδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * Φ_xxαα" == p.to_string() # 1st step:  Φ_xxωω  (here no other greek symbols -> simplified to alpha)
        p.clean_up(set_to_zero=True)
        assert "+ 2 * Φ_xxαα" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "x", "delta", "x"])], prefactor=3)
        p.sort_elements()
        assert "+ 3 * Φ_xxαδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 3 * Φ_xxαα" == p.to_string()# 1st step:  Φ_xxωω  (here no other greek symbols -> simplified to alpha)

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "beta", "delta", "beta"])], prefactor=3)
        p.sort_elements()
        assert "+ 3 * Φ_xββδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 3 * Φ_xxββ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "delta", "beta"])], prefactor=3)
        p.sort_elements()
        assert "+ 3 * Φ_αββδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 3 * Φ_ααββ" == p.to_string()  # Φ_ααββ or Φ_ββδδ

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha", "delta", "gamma"])], prefactor=3)
        p.sort_elements()
        assert "+ 3 * Φ_ααγδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 3 * Φ_ααββ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x", "z", "delta", "beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * Φ_xzβδ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * Φ_xzβδ" == p.to_string()  # β = z, δ = x OR β = x, δ = z -> Factor 2 for term -> no do-able for single Multipole, needs to consider other factors -> not functionality of simplify
        p.clean_up(set_to_zero=True)
        assert "+ 2 * Φ_xzβδ" == p.to_string()


    def test_SimplifyMultipole(self):
        ### MIXED TERMS ###
        p = ProductTerm()
        p.set_elements( [ MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["delta"]) ], prefactor=2)
        p.sort_elements()
        assert "+ 2 * μ_δ * Θ_αβ" == p.to_string()
        p.simplify_Multipole(max_order=1)
        assert "+ 2 * μ_y * Θ_αβ" == p.to_string()
        p.simplify_Multipole(max_order=2)
        assert "+ 2 * μ_y * Θ_αα" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * μ_y * Θ_αα" == p.to_string()

        p = ProductTerm()
        p.set_elements( [ MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["beta"]) ], prefactor=2)
        p.sort_elements()
        assert "+ 2 * μ_β * Θ_αβ" == p.to_string()
        p.simplify_Multipole(max_order=1)
        assert "+ 2 * μ_y * Θ_yα" == p.to_string()
        p.simplify_Multipole(max_order=2)
        assert "+ 2 * μ_y * Θ_yy" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * μ_y * Θ_yy" == p.to_string()

        p = ProductTerm()
        p.set_elements( [ R("delta"), MultipoleMoment(["alpha", "x"]), MultipoleMoment(["delta"]) ], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_δ * μ_δ * Θ_xα" == p.to_string()
        p.simplify_Multipole(max_order=1)
        assert "+ 2 * R_y * μ_y * Θ_xα" == p.to_string()
        p.simplify_Multipole(max_order=2)
        assert "+ 2 * R_y * μ_y * Θ_xx" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * R_y * μ_y * Θ_xx" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("delta"), MultipoleMoment(["x", "x"]), MultipoleMoment(["x", "x", "z","beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_δ * Θ_xx * Φ_xxzβ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * R_δ * Θ_xx * Φ_xxzz" == p.to_string()
        p.clean_up(set_to_zero=True)
        assert "+ 2 * R_δ * Θ_xx * Φ_xxzz" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("delta"), MultipoleMoment(["x", "x"]), MultipoleMoment(["beta", "x", "z", "alpha"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_δ * Θ_xx * Φ_xzαβ" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * R_δ * Θ_xx * Φ_xzαβ" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("delta"), MultipoleMoment(["x", "beta"]), MultipoleMoment(["beta", "x", "z", "alpha"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_δ * Θ_xβ * Φ_xzαβ" == p.to_string()
        p.simplify_Multipole(max_order=2)
        assert "+ 2 * R_δ * Θ_xx * Φ_xxzα" == p.to_string()
        p.simplify_Multipole()
        assert "+ 2 * R_δ * Θ_xx * Φ_xxzz" == p.to_string()



    def test_reduce_indices(self):
        ### BBSE 1: nothing changes
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha"])], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * μ_α" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("alpha")], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * R_α" == p.to_string()

        p = ProductTerm()
        p.set_elements([Delta("alpha", "beta")], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * delta(α, β)" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "gamma"])], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * Ω_αβγ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta"]), Delta("alpha", "beta"), R("beta")], prefactor=2)
        p.clean_up()# mainly sorting elements
        assert "+ 2 * R_β * Θ_αβ * delta(α, β)" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * R_β * Θ_αβ * delta(α, β)" == p.to_string()

        ### BBSE 2: one element changes:
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["beta"])], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * μ_α" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "gamma", "gamma"])], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * Ω_αββ" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("gamma")], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * R_α" == p.to_string()

        p = ProductTerm()
        p.set_elements([Delta("gamma", "beta")], prefactor=2)
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * delta(α, β)" == p.to_string()

        factories = [
            lambda: R("delta"),
            lambda: MultipoleMoment(["zeta", "epsilon"]),
            lambda: Delta("omega", "epsilon"),
        ]
        for order in permutations(factories):
            p = ProductTerm()
            p.set_elements([f() for f in order] , prefactor=2)
            p.sort_elements()# ! no clean up <-> terms might cancel due to system conditions
            assert "+ 2 * R_δ * Θ_εζ * delta(ε, ω)" == p.to_string()
            p.reduce_indices(consider_index_sorting=False)
            assert "+ 2 * R_δ * Θ_αβ * delta(α, γ)" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y", "z", "delta", "delta"])], prefactor=-4)
        p.sort_elements()
        assert "- 4 * Φ_yzδδ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "- 4 * Φ_yzαα" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y", "z", "delta", "delta"]), R("z"), MultipoleMoment(["y"])], prefactor=-4)
        p.sort_elements()
        assert "- 4 * R_z * μ_y * Φ_yzδδ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "- 4 * R_z * μ_y * Φ_yzαα" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y", "y", "delta", "delta"]), R("z"), MultipoleMoment(["y"])], prefactor=-4)
        p.sort_elements()
        assert "- 4 * R_z * μ_y * Φ_yyδδ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "- 4 * R_z * μ_y * Φ_yyαα" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y", "y", "delta", "alpha"]), R("z"), MultipoleMoment(["y"])], prefactor=-4)
        p.sort_elements()
        assert "- 4 * R_z * μ_y * Φ_yyαδ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "- 4 * R_z * μ_y * Φ_yyαβ" == p.to_string()
        p.simplify_Multipole()
        assert "- 4 * R_z * μ_y * Φ_yyαα" == p.to_string()

        # BBSE 3: coordinates included
        p = ProductTerm()
        p.set_elements([R("z"), MultipoleMoment(["y"]), MultipoleMoment(["y", "beta", "beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_z * μ_y * Ω_yββ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * R_z * μ_y * Ω_yαα" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("x"), MultipoleMoment(["y"]), MultipoleMoment(["z", "beta", "gamma"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_x * μ_y * Ω_zβγ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * R_x * μ_y * Ω_zαβ" == p.to_string()

        p = ProductTerm()
        p.set_elements([Delta("x", "epsilon"), MultipoleMoment(["z","gamma"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * Θ_zγ * delta(x, ε)" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * Θ_zα * delta(x, β)" == p.to_string()

        p = ProductTerm()
        p.set_elements([R("z"), MultipoleMoment(["z", "z", "gamma"]), MultipoleMoment(["z", "delta", "delta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * R_z * Ω_zzγ * Ω_zδδ" == p.to_string()
        p.reduce_indices(consider_index_sorting=False)
        assert "+ 2 * R_z * Ω_zzα * Ω_zββ" == p.to_string()


    def test_has_index(self):
        m = MultipoleMoment(["alpha"])
        assert m.has_index(Index("alpha"))
        assert m.has_index(Index("α"))

        r = R("alpha")
        assert r.has_index(Index("alpha"))
        assert r.has_index(Index("α"))

        d = Delta("alpha", "beta")
        assert d.has_index(Index("alpha"))
        assert d.has_index(Index("beta"))
        assert d.has_index(Index("α"))
        assert d.has_index(Index("β"))

    def test_replace_index(self):
        m = MultipoleMoment(["alpha"])
        assert m.to_string() == "μ_α"
        m.replace_index(to_be_replaced = Index("alpha"), replacement = Index("x"))
        assert m.to_string() == "μ_x"

        m = MultipoleMoment(["alpha"])
        assert m.to_string() == "μ_α"
        m.replace_index(to_be_replaced=Index("alpha"), replacement=Index("β"))
        assert m.to_string() == "μ_β"


    def test_reduce_indices_sorting(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["alpha"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * μ_α * Θ_αβ" == p.to_string()
        p.reduce_indices_sorting()
        assert "+ 2 * μ_α * Θ_αβ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["beta", "alpha"]), MultipoleMoment(["beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * μ_β * Θ_αβ" == p.to_string()
        p.reduce_indices_sorting()
        assert "+ 2 * μ_α * Θ_αβ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["z", "z", "beta"]), MultipoleMoment(["alpha", "alpha", "beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * Ω_zzβ * Ω_ααβ" == p.to_string()
        p.reduce_indices_sorting()
        assert "+ 2 * Ω_zzα * Ω_αββ" == p.to_string()

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["z", "z", "alpha"]), MultipoleMoment(["alpha", "beta", "beta"])], prefactor=2)
        p.sort_elements()
        assert "+ 2 * Ω_zzα * Ω_αββ" == p.to_string()
        p.reduce_indices_sorting()
        assert "+ 2 * Ω_zzα * Ω_αββ" == p.to_string()

        p = ProductTerm()
        m1 = MultipoleMoment(["y", "beta", "beta"])
        m1.molecule = "B"
        m2 = MultipoleMoment(["y", "alpha", "alpha"])
        m2.molecule = "A"
        p.set_elements([R("z"), m1, m2 ], prefactor=45)
        assert p.to_string() == "+ 45 * R_z * Ω^B_yββ * Ω^A_yαα"
        p.reduce_indices_sorting()
        assert p.to_string() ==  "+ 45 * R_z * Ω^B_yαα * Ω^A_yββ"


    def test_traceless_conditions(self):
        m = MultipoleMoment(["alpha"])
        traceless = m.traceless_condition()
        assert len(traceless) == 0


        m = MultipoleMoment(["alpha", "beta"])
        assert len(m.traceless_condition()) == 0

        m = MultipoleMoment(["alpha", "alpha"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert traceless[0] == Index("alpha")

        m = MultipoleMoment(["beta", "beta"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert traceless[0] == Index("beta")

        m = MultipoleMoment(["alpha", "beta", "gamma"])
        assert len(m.traceless_condition()) == 0# this sums up all existing components -> 3 * xxy + 1 * yyy + 3 * yzz

        m = MultipoleMoment(["alpha", "beta", "beta"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert Index("beta") in traceless

        m = MultipoleMoment(["alpha", "beta", "z"])
        assert len(m.traceless_condition()) == 0

        m = MultipoleMoment(["beta", "beta", "z"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert Index("beta") in traceless

        m = MultipoleMoment(["y", "alpha", "alpha"])
        assert m.to_string() == f"Ω_yαα"
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert Index("alpha") in traceless

        m = MultipoleMoment(["alpha", "beta", "gamma", "delta"])
        traceless = m.traceless_condition()
        assert len(traceless) == 0

        m = MultipoleMoment(["alpha", "beta", "gamma", "delta"])
        traceless = m.traceless_condition()
        assert len(traceless) == 0

        m = MultipoleMoment(["alpha", "alpha", "delta", "delta"])
        traceless = m.traceless_condition()
        assert len(traceless) == 2
        assert Index("alpha") in traceless
        assert Index("delta") in traceless

        m = MultipoleMoment(["x", "alpha", "x", "alpha"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert Index("alpha") in traceless

        m = MultipoleMoment(["y", "z", "z"])
        traceless = m.traceless_condition()
        assert len(traceless) == 0

        m = MultipoleMoment(["y", "alpha", "alpha"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert Index("alpha") in traceless

        m = MultipoleMoment(["beta", "beta", "x"])
        traceless = m.traceless_condition()
        assert len(traceless) == 1
        assert Index("beta") in traceless

    def test_use_traceless_conditions(self):
        m = MultipoleMoment(["alpha", "alpha"])
        p = ProductTerm()
        p.set_elements([m, R("alpha"), R("beta")], prefactor = 1)
        p.clean_up()
        assert p.to_string() == "+ 1 * R_α * R_β * Θ_αα"
        p.use_traceless_conditions()# traceless condition applies, but cannot be applied due to other coupled terms
        assert p.to_string() == "+ 1 * R_α * R_β * Θ_αα"

        m = MultipoleMoment(["alpha", "alpha"])
        p = ProductTerm()
        p.set_elements([m, R("beta"), R("beta")], prefactor = 1)
        p.clean_up()
        assert p.to_string() == "+ 1 * R_β^(2) * Θ_αα"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        m = MultipoleMoment(["beta", "beta"])
        p = ProductTerm()
        p.set_elements([m, R("alpha"), R("gamma")], prefactor = 1)
        p.clean_up()
        assert p.to_string() == "+ 1 * R_α * R_γ * Θ_ββ"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        m = MultipoleMoment(["beta", "beta"])
        p = ProductTerm()
        p.set_elements([m, MultipoleMoment(["alpha", "gamma"])], prefactor = 1)
        p.clean_up()
        assert p.to_string() == "+ 1 * Θ_αγ * Θ_ββ"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        m = MultipoleMoment(["beta", "beta", "beta"])
        p = ProductTerm()
        p.set_elements([m, MultipoleMoment(["alpha", "gamma"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Ω_βββ * Θ_αγ"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 1 * Ω_βββ * Θ_αγ"


        m = MultipoleMoment(["beta", "beta", "delta"])
        p = ProductTerm()
        p.set_elements([m, MultipoleMoment(["alpha", "delta"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Ω_ββδ * Θ_αδ"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 0"

        m = MultipoleMoment(["beta", "gamma", "gamma"])
        p = ProductTerm()
        p.set_elements([m, MultipoleMoment(["alpha", "gamma"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Ω_βγγ * Θ_αγ"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 1 * Ω_βγγ * Θ_αγ"

        m = MultipoleMoment(["beta", "beta", "x"])
        p = ProductTerm()
        p.set_elements([m, MultipoleMoment(["alpha", "delta"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Ω_xββ * Θ_αδ"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 0"
        for x in ["z", "x"]:
            m = MultipoleMoment(["beta", "beta", x])
            p = ProductTerm()
            p.set_elements([m, MultipoleMoment(["alpha", "delta"])], prefactor=1)
            assert p.to_string() ==  f"+ 1 * Ω_{x}ββ * Θ_αδ"
            p.use_traceless_conditions()
            assert p.to_string() == "+ 0"

        m = MultipoleMoment(["y", "z", "z"])
        p = ProductTerm()
        p.set_elements([m, MultipoleMoment(["y", "alpha", "alpha"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Ω_yzz * Ω_yαα"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 0"

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["x"]), MultipoleMoment(["x", "beta", "beta"]), R("z") ], prefactor=3)
        assert p.to_string() == f"+ 3 * μ_x * Ω_xββ * R_z"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 0"# index y of traceless condition is used in other places, however, it is already a coordinate -> traceless condition still applicable


        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha"]), MultipoleMoment(["alpha"])], prefactor=1)
        assert p.to_string() == f"+ 1 * μ_α * μ_α"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 1 * μ_α * μ_α"

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha","alpha"]), MultipoleMoment(["alpha", "alpha"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Θ_αα * Θ_αα"
        p.use_traceless_conditions()# traceless condition would apply for separate terms, but cannot be used since sum is coupled
        assert p.to_string() == f"+ 1 * Θ_αα * Θ_αα"

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["y"]), MultipoleMoment(["y", "alpha", "alpha"])], prefactor=1)
        assert p.to_string() == f"+ 1 * μ_y * Ω_yαα"
        p.use_traceless_conditions()
        assert p.to_string() == f"+ 0"


        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha","alpha"]), MultipoleMoment(["alpha","alpha", "x", "x"])], prefactor=1)
        p.use_traceless_conditions()
        assert p.to_string() == "+ 1 * Θ_αα * Φ_xxαα"

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "alpha", "beta", "beta"])], prefactor=1)
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"  # because beta itself is sufficient, even though the alpha sum does not necessarily vanish



    def test_find_unresolved_multipoles(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "alpha"])], prefactor=1)
        assert p.to_string() == f"+ 1 * Θ_αα * Θ_αα"
        result = p.includes_duplicates()
        assert result[0] is True
        assert len(result[1]) == 1
        assert result[1][0] == Index("alpha")

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha"]), MultipoleMoment(["alpha"])], prefactor=1)
        result = p.includes_duplicates()
        assert result[0] is True
        assert len(result[1]) == 1
        assert result[1][0] == Index("alpha")


        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "beta"])], prefactor=1)
        result = p.includes_duplicates()
        assert result[0] is False
        assert len(result[1]) == 0

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "alpha", "alpha"])], prefactor=1)
        result = p.includes_duplicates()
        assert result[0] is False
        assert len(result[1]) == 0








