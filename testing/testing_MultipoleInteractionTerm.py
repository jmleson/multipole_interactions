import unittest
import sympy as sp

from SummandTerm import MultipoleInteraction


class TestMultipoleInteractionTerm(unittest.TestCase):


    def test_set_r_prefactor(self):
        n = sp.symbols('n')

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        assert s.order == 2
        assert s.prefactor_expansion == 1/n.subs(n, 1)
        assert s.r_prefactor.exponent == -5

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        assert s.order == 3
        assert s.prefactor_expansion == -1 / n.subs(n, 3)
        assert s.r_prefactor.exponent == -7

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        assert s.order == 4
        assert s.prefactor_expansion == 1/n.subs(n, 15)
        assert s.r_prefactor.exponent == -9

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=4)
        assert s.order == 5
        assert s.prefactor_expansion == -1/n.subs(n, 105)
        assert s.r_prefactor.exponent == -11


        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=2)
        assert s.order == 4
        assert s.prefactor_expansion == 1/n.subs(n, 9)
        assert s.r_prefactor.exponent == -9


    def test_get_name(self):
        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        assert s.get_name() == "T_αβ μ^A_α μ^B_β"

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        assert s.get_name() == "T_αβγ μ^A_α Θ^B_βγ"

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        assert s.get_name() == "T_αβγ Θ^A_αβ μ^B_γ"

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        assert s.get_name() == "T_αβγδ μ^A_α Ω^B_βγδ"

    def test_total(self):
        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "T_αβ μ^A_α μ^B_β = [ - 1 * R_z^(-3) * μ^A_y * μ^B_y ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        s.simplify_in_latex_steps()
        assert "T_αβγ μ^A_α Θ^B_βγ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "T_αβγ Θ^A_αβ μ^B_γ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        s.simplify_in_latex_steps()
        assert "T_αβγδ μ^A_α Ω^B_βγδ = [ - 3 * R_z^(-5) * μ^A_y * Ω^B_yzz ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=4)
        s.simplify_in_latex_steps()
        assert "T_αβγδε μ^A_α Φ^B_βγδε = + 0" == s.full_string_tensor()


        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=2)
        s.simplify_in_latex_steps()
        assert ("T_αβγδ Θ^A_αβ Θ^B_γδ = [ "
                "+ 17/3 * R_z^(-5) * Θ^A_zz * Θ^B_zz "
                "+ 2/3 * R_z^(-5) * Θ^A_xx * Θ^B_xx "
                "+ 2/3 * R_z^(-5) * Θ^A_yy * Θ^B_yy ]") == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=3)
        s.simplify_in_latex_steps()
        assert "T_αβγδε Θ^A_αβ Ω^B_γδε = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=4)
        s.simplify_in_latex_steps()
        assert "T_αβγδεζη Ω^A_αβγ Φ^B_δεζη = + 0" == s.full_string_tensor()



