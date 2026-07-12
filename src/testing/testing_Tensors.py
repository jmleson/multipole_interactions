import unittest

from src.Tensor import Tensor


class TestTensors(unittest.TestCase):
    t0 = Tensor(0)
    t1 = Tensor(1)
    t2 = Tensor(2)
    t3 = Tensor(3)
    t4 = Tensor(4)
    t5 = Tensor(5)
    t6 = Tensor(6)

    def test_set_r_prefactor(self):
        assert self.t0.order == 0
        assert self.t0.r_prefactor.exponent == -1

        assert self.t1.order == 1
        assert self.t1.r_prefactor.exponent == -3

        assert self.t2.order == 2
        assert self.t2.r_prefactor.exponent == -5

        assert self.t3.order == 3
        assert self.t3.r_prefactor.exponent == -7

        assert self.t4.order == 4
        assert self.t4.r_prefactor.exponent == -9

        assert self.t5.order == 5
        assert self.t5.r_prefactor.exponent == -11

        assert self.t6.order == 6
        assert self.t6.r_prefactor.exponent == -13


    def test_get_name(self):
        assert self.t0.get_name() == "T"
        assert self.t1.get_name() == "T_α"
        assert self.t2.get_name() == "T_αβ"
        assert self.t3.get_name() == "T_αβγ"
        assert self.t4.get_name() == "T_αβγδ"

    def test_amount_of_terms(self):
        assert len(self.t1.terms) == 1
        assert len(self.t2.terms) == 2
        assert len(self.t3.terms) == 4
        assert len(self.t4.terms) == 10


    def test_full_string_tensor_basics(self):
        result = self.t0.full_string_tensor_basics()
        assert result[result.find("="):] == "= R^(-1)"

        result = self.t1.full_string_tensor_basics()
        assert result[result.find("="):] == "= R^(-3) * [ - 1 * R_α ]"

        result = self.t2.full_string_tensor_basics()
        assert result[result.find("="):] == "= R^(-5) * [ + 3 * R_α * R_β - 1 * R^(2) * delta(α, β) ]"

        result = self.t3.full_string_tensor_basics()
        assert result[result.find("="):] == ("= R^(-7) * [ - 15 * R_α * R_β * R_γ "
                                             "+ 3 * R^(2) * R_α * delta(β, γ) "
                                             "+ 3 * R^(2) * R_β * delta(α, γ) "
                                             "+ 3 * R^(2) * R_γ * delta(α, β) ]")

        result = self.t4.full_string_tensor_basics()
        assert result[result.find("="):] == ("= R^(-9) * [ "
                                             "+ 105 * R_α * R_β * R_γ * R_δ "
                                             "- 15 * R^(2) * R_α * R_β * delta(γ, δ) "
                                             "- 15 * R^(2) * R_α * R_γ * delta(β, δ) "
                                             "- 15 * R^(2) * R_α * R_δ * delta(β, γ) "
                                             "- 15 * R^(2) * R_β * R_γ * delta(α, δ) "
                                             "- 15 * R^(2) * R_β * R_δ * delta(α, γ) "
                                             "- 15 * R^(2) * R_γ * R_δ * delta(α, β) "
                                             "+ 3 * R^(4) * delta(α, β) * delta(γ, δ) "
                                             "+ 3 * R^(4) * delta(α, γ) * delta(β, δ) "
                                             "+ 3 * R^(4) * delta(α, δ) * delta(β, γ) ]")


