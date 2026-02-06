import unittest

from SummandTerm import MultipoleInteraction


class TestMultipoleInteractionTerm(unittest.TestCase):


    def test_set_r_prefactor(self):
        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=1)
        assert s.order == 2
        assert s.r_prefactor.exponent == -5

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        assert s.order == 3
        assert s.r_prefactor.exponent == -7

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=3)
        assert s.order == 4
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
        s.simplify()
        assert "T_αβ μ_α μ_β =  [ + 1 * R_z^(-3) * μ_y * μ_y ]" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=1, multipole_order_2=2)
        s.simplify()
        assert "T_αβγ μ_α Θ_βγ = + 0" == s.full_string_tensor()

        s = MultipoleInteraction(multipole_order_1=2, multipole_order_2=1)
        s.simplify()
        assert "T_αβγ Θ_αβ μ_γ = + 0" == s.full_string_tensor()


