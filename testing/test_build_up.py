import unittest
import sympy as sp

from Multipole import MultipoleMoment
from build_up_tensor_equations.construct_full_tensor import r_zero_part, r_two_part
from build_up_tensor_equations.tensor import TensorTerm, Tensor
from build_up_tensor_equations.variations_between_R_and_delta import variations_between_R_and_delta


class TestBuildUp(unittest.TestCase):

    def test_Tensor(self):
        dipole = MultipoleMoment(indices=["alpha"])
        assert dipole.prefactor_in_potential == -1
        t = Tensor(dipole)
        assert t.order == 1
        assert t.get_prefactor() == sp.sympify("-1")
        assert t.get_tensor_notation() == "T_α"

        quadrupole = MultipoleMoment(indices=["alpha", "beta"])
        assert quadrupole.prefactor_in_potential == 1
        t = Tensor(quadrupole)
        assert t.order == 2
        assert t.get_prefactor() == sp.sympify("1/3")
        assert t.get_tensor_notation() == "T_αβ"

        oktopole = MultipoleMoment(indices=["alpha", "beta", "gamma"])
        assert oktopole.prefactor_in_potential == -1
        t = Tensor(oktopole)
        assert t.order == 3
        assert t.get_prefactor() == sp.sympify("-1/15")
        assert t.get_tensor_notation() == "T_αβγ"


    def test_tensor_to_equation(self):
        dipole1 = MultipoleMoment(indices=["alpha"])
        dipole2 = MultipoleMoment(indices=["beta"])
        t = TensorTerm(m1=dipole1, m2=dipole2)
        assert t.order == 2
        assert t.get_prefactor() == sp.sympify("-1")
        assert "-1 * T_αβ * μ_α * μ_β" == t.to_string()

        dipole = MultipoleMoment(indices=["alpha"])
        quadrupole = MultipoleMoment(indices=["beta", "gamma"])
        t = TensorTerm(m1=dipole, m2=quadrupole)
        assert t.order == 3
        assert t.get_prefactor() == sp.sympify("1/3")
        assert "1/3 * T_αβγ * μ_α * Θ_βγ" == t.to_string()

        dipole = MultipoleMoment(indices=["alpha"])
        quadrupole = MultipoleMoment(indices=["beta", "gamma", "delta"])
        t = TensorTerm(m1=dipole, m2=quadrupole)
        assert t.order == 4
        assert t.get_prefactor() == sp.sympify("-1/15")
        assert "-1/15 * T_αβγδ * μ_α * Ω_βγδ" == t.to_string()

        dipole = MultipoleMoment(indices=["alpha", "beta"])
        quadrupole = MultipoleMoment(indices=["gamma", "delta"])
        t = TensorTerm(m1=dipole, m2=quadrupole)
        assert t.order == 4
        assert t.get_prefactor() == sp.sympify("1/9")
        assert "1/9 * T_αβγδ * Θ_αβ * Θ_γδ" == t.to_string()


    def test_r_zero_part(self):
        part = r_zero_part(order=1)
        assert "-1 * R_α" == part.to_string()

        part = r_zero_part(order=2)
        assert "3 * R_α * R_β" == part.to_string()

        part = r_zero_part(order=3)
        assert "-15 * R_α * R_β * R_γ" == part.to_string()

        part = r_zero_part(order=4)
        assert "105 * R_α * R_β * R_γ * R_δ" == part.to_string()

        part = r_zero_part(order=5)
        assert "-945 * R_α * R_β * R_γ * R_δ * R_ε" == part.to_string()

    def test_r_two_part(self):
        part = r_two_part(order=2)
        assert "-1 * R_z^(2) * delta(α, β)" == part.to_string()

        part = r_two_part(order=3)
        assert "3 * R_z^(2) * R_α * delta(β, γ)" == part.to_string()

        part = r_two_part(order=4)
        assert "-15 * R_z^(2) * R_α * R_β * delta(γ, δ)" == part.to_string()

    def test_variations_between_R_and_delta(self):
        result = variations_between_R_and_delta(1,1)
        assert len(result) == 3
        strings = [p.to_string() for p in result]
        assert "1 * R_α * delta(β, γ)" in strings
        assert "1 * R_β * delta(α, γ)" in strings
        assert "1 * R_γ * delta(α, β)" in strings

        result = variations_between_R_and_delta(2, 1)
        assert len(result) == 6
        strings = [p.to_string() for p in result]
        assert "1 * R_α * R_β * delta(γ, δ)" in strings
        assert "1 * R_γ * R_δ * delta(α, β)" in strings
        assert "1 * R_β * R_γ * delta(α, δ)" in strings
        assert "1 * R_β * R_δ * delta(α, γ)" in strings
        assert "1 * R_α * R_γ * delta(β, δ)" in strings
        assert "1 * R_α * R_δ * delta(β, γ)" in strings

