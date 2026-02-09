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
        assert s.get_name() == "T_αβ μ_α μ_β"

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        assert s.get_name() == "T_αβγ μ_α Θ_βγ"

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        assert s.get_name() == "T_αβγ Θ_αβ μ_γ"

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        assert s.get_name() == "T_αβγδ μ_α Ω_βγδ"

    def test_total(self):
        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "T_αβ μ_α μ_β = [ - 1 * R_z^(-3) * μ_y * μ_y ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        s.simplify_in_latex_steps()
        assert "T_αβγ μ_α Θ_βγ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        s.simplify_in_latex_steps()
        assert "T_αβγ Θ_αβ μ_γ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        s.simplify_in_latex_steps()
        assert "T_αβγδ μ_α Ω_βγδ = [ - 3 * R_z^(-5) * μ_y * Ω_yzz ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=4)
        s.simplify_in_latex_steps()
        assert "T_αβγδε μ_α Φ_βγδε = + 0" == s.full_string_tensor()


        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=2)
        s.simplify_in_latex_steps()
        assert ("T_αβγδ Θ_αβ Θ_γδ = [ "
                "+ 17/3 * R_z^(-5) * Θ_zz * Θ_zz "
                "+ 2/3 * R_z^(-5) * Θ_xx * Θ_xx "
                "+ 2/3 * R_z^(-5) * Θ_yy * Θ_yy ]") == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=3)
        s.simplify_in_latex_steps()
        assert "T_αβγδε Θ_αβ Ω_γδε = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=3, multipole_order_2=4)
        s.simplify_in_latex_steps()
        assert "T_αβγδεζη Ω_αβγ Φ_δεζη = + 0" == s.full_string_tensor()



