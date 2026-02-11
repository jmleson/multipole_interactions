import unittest
import sympy as sp

from SummandTerm import MultipoleInteraction
from term_components.Multipole import MultipoleMoment
from term_components.ProductTerm import ProductTerm
from term_components.get_product_terms_from_greek_sum import get_product_terms_from_greek_sum_all


class TestMultipoleInteractionTerm(unittest.TestCase):


    def test_set_r_prefactor(self):
        n = sp.symbols('n')

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        assert s.order == 2
        assert s.prefactor_of_expansion == - 1 / n.subs(n, 1)
        assert s.r_prefactor.exponent == -5

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        assert s.order == 3
        assert s.prefactor_of_expansion == +1 / n.subs(n, 3)
        assert s.r_prefactor.exponent == -7

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        assert s.order == 3
        assert s.prefactor_of_expansion == -1 / n.subs(n, 3)
        assert s.r_prefactor.exponent == -7

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        assert s.order == 4
        assert s.prefactor_of_expansion == -1 / n.subs(n, 15)
        assert s.r_prefactor.exponent == -9

        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=1)
        assert s.order == 4
        assert s.prefactor_of_expansion == -1 / n.subs(n, 15)# ! here identical to switched version
        assert s.r_prefactor.exponent == -9

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=4)
        assert s.order == 5
        assert s.prefactor_of_expansion == +1 / n.subs(n, 105)
        assert s.r_prefactor.exponent == -11


        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=2)
        assert s.order == 4
        assert s.prefactor_of_expansion == +1 / n.subs(n, 9)
        assert s.r_prefactor.exponent == -9

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=3)
        assert s.order == 5
        assert s.prefactor_of_expansion == -1 / n.subs(n, 3*15)
        assert s.r_prefactor.exponent == -11

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=4)
        assert s.order == 6
        assert s.prefactor_of_expansion == 1 / n.subs(n, 3 * 105)
        assert s.r_prefactor.exponent == -13

        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=3)
        assert s.order == 6
        assert s.prefactor_of_expansion == -1 / n.subs(n, 15*15)
        assert s.r_prefactor.exponent == -13


    def test_get_name(self):
        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        assert s.get_name() == "- 1 * T_αβ μ^B_α μ^A_β"

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        assert s.get_name() == "+ 1/3 * T_αβγ μ^B_α Θ^A_βγ"

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        assert s.get_name() == "- 1/3 * T_αβγ Θ^B_αβ μ^A_γ"

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        assert s.get_name() == "- 1/15 * T_αβγδ μ^B_α Ω^A_βγδ"

        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=1)
        assert s.get_name() == "- 1/15 * T_αβγδ Ω^B_αβγ μ^A_δ"

    def test_find_unresolved_multipoles(self):
        m = MultipoleInteraction(multipole_order_1 = 1, multipole_order_2 = 1)
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha"]), MultipoleMoment(["alpha"])], prefactor=1)
        m.tensor_terms = [p]
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= R_z^(-5) * -1 * [ + 1 * μ_α * μ_α ]"
        m.find_unresolved_multipoles(get_product_terms_from_greek_sum_all)
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= R_z^(-5) * -1 * [ + 1 * μ_x * μ_x + 1 * μ_y * μ_y + 1 * μ_z * μ_z ]"


        m = MultipoleInteraction(multipole_order_1 = 3, multipole_order_2 = 3)
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "x", "x"])], prefactor=1)
        m.tensor_terms = [p]
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= R_z^(-13) * -1/225 * [ + 1 * Ω_xxα ]"
        m.find_unresolved_multipoles(get_product_terms_from_greek_sum_all)
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= R_z^(-13) * -1/225 * [ + 1 * Ω_xxx + 1 * Ω_xxy + 1 * Ω_xxz ]"
        m.clean_up()
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= [ - 1/225 * R_z^(-13) * Ω_xxy ]"


        m = MultipoleInteraction(multipole_order_1 = 3, multipole_order_2 = 3)
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "beta"])], prefactor=1)
        m.tensor_terms = [p]
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= R_z^(-13) * -1/225 * [ + 1 * Ω_αββ ]"
        m.find_unresolved_multipoles(get_product_terms_from_greek_sum_all)
        result = m.full_string_tensor()
        assert (result[result.find("="):] == ("= R_z^(-13) * -1/225 * [ "
                                             "+ 1 * Ω_xxx + 1 * Ω_xxy + 1 * Ω_xxz "
                                             "+ 1 * Ω_xyy + 1 * Ω_yyy + 1 * Ω_yyz "
                                             "+ 1 * Ω_xzz + 1 * Ω_yzz + 1 * Ω_zzz ]")
            or result[result.find("="):] == ("= R_z^(-13) * -1/225 * [ "
                                             "+ 1 * Ω_xxx + 1 * Ω_xyy + 1 * Ω_xzz "
                                             "+ 1 * Ω_xxy + 1 * Ω_yyy + 1 * Ω_yzz "
                                             "+ 1 * Ω_xxz + 1 * Ω_yyz + 1 * Ω_zzz ]")
                )
        m.clean_up()
        result = m.full_string_tensor()
        assert result[result.find("="):] == "= [ - 1/225 * R_z^(-13) * Ω_xxy - 1/225 * R_z^(-13) * Ω_yyy - 1/225 * R_z^(-13) * Ω_yzz ]"



    def test_total(self):
        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "- 1 * T_αβ μ^B_α μ^A_β = [ + 1 * R_z^(-3) * μ^B_y * μ^A_y ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        s.simplify_in_latex_steps()
        assert "+ 1/3 * T_αβγ μ^B_α Θ^A_βγ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "- 1/3 * T_αβγ Θ^B_αβ μ^A_γ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        s.simplify_in_latex_steps()
        assert "- 1/15 * T_αβγδ μ^B_α Ω^A_βγδ = [ + 3 * R_z^(-5) * μ^B_y * Ω^A_yzz ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "- 1/15 * T_αβγδ Ω^B_αβγ μ^A_δ = [ + 3 * R_z^(-5) * Ω^B_yzz * μ^A_y ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=4)
        s.simplify_in_latex_steps()
        assert "+ 1/105 * T_αβγδε μ^B_α Φ^A_βγδε = + 0" == s.full_string_tensor()



        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=2)
        s.simplify_in_latex_steps()
        assert ("+ 1/9 * T_αβγδ Θ^B_αβ Θ^A_γδ = [ "
                "+ 17/3 * R_z^(-5) * Θ^B_zz * Θ^A_zz "
                "+ 2/3 * R_z^(-5) * Θ^B_xx * Θ^A_xx "
                "+ 2/3 * R_z^(-5) * Θ^B_yy * Θ^A_yy ]") == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=3)
        s.simplify_in_latex_steps()
        assert "- 1/45 * T_αβγδε Θ^B_αβ Ω^A_γδε = + 0" == s.full_string_tensor()


        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=4)
        s.simplify_in_latex_steps()
        assert "+ 1/1575 * T_αβγδεζη Ω^B_αβγ Φ^A_δεζη = + 0" == s.full_string_tensor()



