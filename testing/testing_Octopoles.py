import copy
import unittest

from term_components.Multipole import MultipoleMoment
from term_components.ProductTerm import ProductTerm
from term_components.R import R


class testingOctopoles(unittest.TestCase):
    # 0 identical:
    m1 = MultipoleMoment(["alpha", "beta", "gamma"])# = 3 * xxy + 1 * yyy + 3 * yzz
    # 2 identical:
    m2 = MultipoleMoment(["alpha", "beta", "beta"])# = y beta beta = 0
    # 3 identical:
    m3 = MultipoleMoment(["alpha", "alpha", "alpha"])# = yyy

    m4 = MultipoleMoment(["x", "alpha", "alpha"])# always 0
    m5 = MultipoleMoment(["x", "alpha", "beta"])# = 2 * xxy
    m6 = MultipoleMoment(["y", "alpha", "alpha"])# = xxy + yyy + yzz -> 0
    m7 = MultipoleMoment(["y", "alpha", "beta"])# = xxy + yyy + yzz -> 0

    m8 = MultipoleMoment(["x", "x", "alpha"])# = xxy
    m9 = MultipoleMoment(["x", "y", "alpha"])# = xxy
    m10 = MultipoleMoment(["x", "z", "alpha"])# always 0
    m11 = MultipoleMoment(["y", "y", "alpha"])# = yyy

    o1 = MultipoleMoment(["x", "x", "x"])
    o2 = MultipoleMoment(["x", "x", "y"])
    o3 = MultipoleMoment(["x", "x", "z"])
    o4 = MultipoleMoment(["x", "y", "y"])
    o5 = MultipoleMoment(["x", "z", "z"])
    o6 = MultipoleMoment(["y", "y", "y"])
    o7 = MultipoleMoment(["y", "y", "z"])
    o8 = MultipoleMoment(["y", "z", "z"])
    o9 = MultipoleMoment(["z", "z", "z"])
    o10 = MultipoleMoment(["x", "y", "z"])


    def test_set_up(self):#
        def get_simple_copy(m):
            p = ProductTerm()
            p.set_elements([copy.deepcopy(m)], prefactor=1)
            p.simplify_Multipole()
            p.use_traceless_conditions()
            p.clean_up()
            return p

        assert get_simple_copy(self.m2).to_string() == "+ 0"
        assert get_simple_copy(self.m3).to_string() == "+ 1 * " + self.o6.to_string()
        assert get_simple_copy(self.m4).to_string() == "+ 0"
        assert get_simple_copy(self.m6).to_string() == "+ 0"
        assert get_simple_copy(self.m7).to_string() == "+ 0"
        assert get_simple_copy(self.m8).to_string() == "+ 1 * " + self.o2.to_string()
        assert get_simple_copy(self.m9).to_string() == "+ 1 * " + self.o2.to_string()
        assert get_simple_copy(self.m10).to_string() == "+ 0"
        assert get_simple_copy(self.m11).to_string() == "+ 1 * " + self.o6.to_string()


    def test_is_zero_for_whatever_reason(self):
        # zero-components:
        assert self.m2.is_zero() or len(self.m2.traceless_condition()) > 0
        assert self.m4.is_zero() and len(self.m4.traceless_condition()) > 0
        assert self.m6.is_zero() or len(self.m6.traceless_condition()) > 0
        assert self.m10.is_zero() or len(self.m10.traceless_condition()) > 0
        assert self.o1.is_zero()
        assert self.o3.is_zero()
        assert self.o4.is_zero()
        assert self.o5.is_zero()
        assert self.o7.is_zero()
        assert self.o9.is_zero()
        assert self.o10.is_zero()

        # non-zero components:
        assert not self.m1.is_zero() and len(self.m1.traceless_condition()) == 0
        assert not self.m3.is_zero() and len(self.m3.traceless_condition()) == 0
        assert not self.m5.is_zero() and len(self.m5.traceless_condition()) == 0
        assert not self.m7.is_zero() and len(self.m7.traceless_condition()) == 0
        assert not self.m8.is_zero() and len(self.m8.traceless_condition()) == 0
        assert not self.m9.is_zero() and len(self.m9.traceless_condition()) == 0
        assert not self.m11.is_zero() and len(self.m11.traceless_condition()) == 0
        assert not self.o2.is_zero() and len(self.o2.traceless_condition()) == 0
        assert not self.o6.is_zero() and len(self.o6.traceless_condition()) == 0
        assert not self.o8.is_zero() and len(self.o8.traceless_condition()) == 0


    def test_simplify(self):
        def get_simple_copy(m):
            p = ProductTerm()
            p.set_elements([copy.deepcopy(m)], prefactor=1)
            p.simplify_Multipole()
            return p.factors[0]

        # unchanged:
        assert self.m1.to_string() == get_simple_copy(self.m1).to_string()# Ω_αβγ -> 2 identical, but different factors
        assert self.m4.to_string() == get_simple_copy(self.m4).to_string()# but is zero
        assert self.m5.to_string() == get_simple_copy(self.m5).to_string()# Ω_xαβ (has to be xxy, but 2x)
        assert self.m6.to_string() == get_simple_copy(self.m6).to_string()# but is zero
        assert self.m10.to_string() == get_simple_copy(self.m10).to_string()# but is zero
        assert self.o1.to_string() == get_simple_copy(self.o1).to_string()
        assert self.o2.to_string() == get_simple_copy(self.o2).to_string()
        assert self.o3.to_string() == get_simple_copy(self.o3).to_string()
        assert self.o4.to_string() == get_simple_copy(self.o4).to_string()
        assert self.o5.to_string() == get_simple_copy(self.o5).to_string()
        assert self.o6.to_string() == get_simple_copy(self.o6).to_string()
        assert self.o7.to_string() == get_simple_copy(self.o7).to_string()
        assert self.o8.to_string() == get_simple_copy(self.o8).to_string()
        assert self.o9.to_string() == get_simple_copy(self.o9).to_string()
        assert self.o10.to_string() == get_simple_copy(self.o10).to_string()

        # changed:
        assert get_simple_copy(self.m2).to_string() == "Ω_yββ"# Ω_αββ
        assert get_simple_copy(self.m3).to_string() == "Ω_yyy"# Ω_ααα -> yyy

        assert get_simple_copy(self.m7).to_string() == "Ω_yαα"# -> y alpha alpha -> dann traceless
        assert len(get_simple_copy(self.m7).traceless_condition()) > 0

        assert get_simple_copy(self.m8).to_string() == "Ω_xxy"# -> xxy
        assert get_simple_copy(self.m9).to_string() == "Ω_xxy"# -> xxy
        assert get_simple_copy(self.m11).to_string() == "Ω_yyy"# -> yyy



    def test_simplifyR(self):
        def get_p(r_index, m):
            r = R(r_index)
            p = ProductTerm()
            p.set_elements([r, copy.deepcopy(m)], prefactor = -2)
            p.sort_elements()
            return p

        p = get_p("alpha", self.m1)
        assert p.to_string() == "- 2 * R_α * Ω_αβγ"
        p.simplify_Multipole()
        assert p.to_string() == "- 2 * R_α * Ω_αβγ"
        p.simplify_R()
        assert p.to_string() == "- 2 * R_z * Ω_zβγ"
        p.simplify_Multipole()
        assert p.to_string() == "- 2 * R_z * Ω_zβγ"# zzy and zyz




    def test_products(self):
        def get_p(m1, m2):
            # self.m2.molecule = "A"
            # self.m1.molecule = "B"
            m1, m2 = copy.deepcopy(m1), copy.deepcopy(m2)
            p = ProductTerm()
            p.set_elements([m1, m2], prefactor=-4)
            p.clean_up()
            p.sort_elements()
            return p


        # combinations of two changeable components:
        p = get_p(self.m2, self.m3)
        assert p.to_string() == "- 4 * Ω_ααα * Ω_αββ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yββ"
        p.reduce_indices(consider_index_sorting=True)
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yαα"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        p = get_p(self.m2, self.m7)
        assert p.to_string() == "- 4 * Ω_yαβ * Ω_αββ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yyy"# "- 4 * Ω_yαβ * Ω_αββ" -> "- 4 * Ω_yyβ * Ω_yββ" -> "- 4 * Ω_yyy * Ω_yyy"

        p = get_p(self.m2, self.m8)
        assert p.to_string() == "- 4 * Ω_xxα * Ω_αββ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_xxy * Ω_yββ"# "- 4 * Ω_xxα * Ω_αββ" -> "- 4 * Ω_xxy * Ω_yββ"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        p = get_p(self.m2, self.m9)
        assert p.to_string() == "- 4 * Ω_xyα * Ω_αββ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_xxy * Ω_xββ"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        p = get_p(self.m2, self.m11)
        assert p.to_string() == "- 4 * Ω_yyα * Ω_αββ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yββ"
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"

        p = get_p(self.m3, self.m7)
        assert p.to_string() == "- 4 * Ω_yαβ * Ω_ααα"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yyy"
        p.use_traceless_conditions()
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yyy"

        # combinations of itself unchangeable components:
        p = get_p(self.m1, self.m1)
        p.simplify_Multipole()
        assert p.to_string() == get_p(self.m1, self.m1).to_string()

        p = get_p(self.m1, self.m4)
        p.simplify_Multipole()
        assert p.to_string() == get_p(self.m1, self.m4).to_string()

        p = get_p(self.m1, self.m5)
        p.simplify_Multipole()
        assert p.to_string() == get_p(self.m1, self.m5).to_string()

        p = get_p(self.m10, self.m6)
        p.simplify_Multipole()
        assert p.to_string() == get_p(self.m10, self.m6).to_string()



        # combinations of one changeable and one not-changeable component:
        p = get_p(self.m1, self.m2)
        assert p.to_string() == "- 4 * Ω_αββ * Ω_αβγ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yαα * Ω_yαα"#  "- 4 * Ω_αββ * Ω_αβγ" -> "- 4 * Ω_yββ * Ω_yβγ" -> "- 4 * Ω_y omega omega * Ω_y omega omega"
        p.use_traceless_conditions()
        assert p.to_string() == "- 4 * Ω_yαα * Ω_yαα"# alpha coupled -> traceless condition not applicable

        p = get_p(self.m1, self.m3)
        assert p.to_string() == "- 4 * Ω_ααα * Ω_αβγ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yyy * Ω_yαα"  # -> "- 4 * Ω_yyy * Ω_yβγ" ->  "- 4 * Ω_yyy * Ω_y omega omega" ->
        p.use_traceless_conditions()
        assert p.to_string() == "+ 0"  #

        p = get_p(self.m7, self.m1)
        assert p.to_string() == "- 4 * Ω_yαβ * Ω_αβγ"
        p.simplify_Multipole()
        assert p.to_string() == "- 4 * Ω_yαα * Ω_yαα"
        p.use_traceless_conditions()
        assert p.to_string() == "- 4 * Ω_yαα * Ω_yαα"

        p = get_p(self.m9, self.m1)
        assert p.to_string() == "- 4 * Ω_xyα * Ω_αβγ"
        p.simplify_Multipole()
        p.reduce_indices(consider_index_sorting=True)
        assert p.to_string() == "- 4 * Ω_xxy * Ω_xαβ"
        p.use_traceless_conditions()
        assert p.to_string() == "- 4 * Ω_xxy * Ω_xαβ"

        for m in [self.m8, self.m11]:
            p = get_p(m, self.m1)
            p.simplify_Multipole()
            p.use_traceless_conditions()
            assert p.to_string() == "+ 0"

        # special case:
        for m in [self.m1, self.m2, self.m3, self.m4, self.m5, self.m6, self.m7, self.m8, self.m9, self.m10, self.m11, self.o1, self.o2, self.o3]:
            p = get_p(m, self.m4)
            assert p.to_string() == "+ 0"# because Ω_xαα is always zero (not only because traceless condition)

            p = get_p(m, self.m10)
            assert p.to_string() == "+ 0"  # because Ω_xαα is always zero (not only because traceless condition)






