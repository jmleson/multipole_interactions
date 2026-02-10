import unittest
from itertools import permutations

import sympy as sp

from term_components.Delta import Delta
from term_components.Index import Index
from term_components.Multipole import MultipoleMoment
from build_up_tensor_equations.construct_full_tensor import construct_all_tensor_terms, r_part
from build_up_tensor_equations.tensor import TensorTerm, Tensor
from build_up_tensor_equations.variations_between_R_and_delta import variations_between_R_and_delta
from term_components.ProductTerm import ProductTerm
from term_components.R import R
from term_components.get_product_terms_from_greek_sum import get_product_terms_from_greek_sum_duplicateMultipoles, \
    get_product_terms_from_greek_sum_all


class TestBuildUp(unittest.TestCase):

    def test_sorting_Index(self):
        i1 = Index("alpha")
        i2 = Index("beta")
        i3 = Index("gamma")
        i4 = Index("delta")

        should_be = [i1,i2,i3,i4]
        for order in list(permutations(should_be)):
            sorted_list = sorted(order)
            assert [i.to_string() for i in sorted_list] == [i.to_string() for i in should_be]


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

        hexadecapole = MultipoleMoment(indices=["alpha", "beta", "gamma", "delta"])
        assert hexadecapole.prefactor_in_potential == 1
        h = Tensor(hexadecapole)
        assert h.order == 4
        assert h.get_prefactor() == sp.sympify("+1/105")
        assert h.get_tensor_notation() == "T_αβγδ"


    def test_tensor_to_equation(self):
        dipole1 = MultipoleMoment(indices=["alpha"])
        dipole2 = MultipoleMoment(indices=["beta"])
        t = TensorTerm(m1=dipole1, m2=dipole2)
        assert t.order == 2
        assert t.get_prefactor() == sp.sympify("-1")
        assert "- 1 * T_αβ * μ_α * μ_β" == t.to_string()

        dipole = MultipoleMoment(indices=["alpha"])
        quadrupole = MultipoleMoment(indices=["beta", "gamma"])
        t = TensorTerm(m1=dipole, m2=quadrupole)
        assert t.order == 3
        assert t.get_prefactor() == sp.sympify("1/3")
        assert "+ 1/3 * T_αβγ * μ_α * Θ_βγ" == t.to_string()

        dipole = MultipoleMoment(indices=["alpha"])
        quadrupole = MultipoleMoment(indices=["beta", "gamma", "delta"])
        t = TensorTerm(m1=dipole, m2=quadrupole)
        assert t.order == 4
        assert t.get_prefactor() == sp.sympify("-1/15")
        assert "- 1/15 * T_αβγδ * μ_α * Ω_βγδ" == t.to_string()

        dipole = MultipoleMoment(indices=["alpha", "beta"])
        dipole.molecule = "A"
        quadrupole = MultipoleMoment(indices=["gamma", "delta"])
        quadrupole.molecule = "B"
        t = TensorTerm(m1=dipole, m2=quadrupole)
        assert t.order == 4
        assert t.get_prefactor() == sp.sympify("1/9")
        assert "+ 1/9 * T_αβγδ * Θ^A_αβ * Θ^B_γδ" == t.to_string()


        dipole = MultipoleMoment(indices=["alpha"])
        hexadecapole = MultipoleMoment(indices=["beta", "gamma", "delta", "epsilon"])
        t = TensorTerm(m1=dipole, m2=hexadecapole)
        assert t.order == 5
        assert t.get_prefactor() == sp.sympify("1/105")
        assert "+ 1/105 * T_αβγδε * μ_α * Φ_βγδε" == t.to_string()

    def test_r_zero_part(self):
        part = r_part(order=0, exponent_of_r=0)
        assert len(part) == 0

        part = r_part(order=1, exponent_of_r=0)
        assert len(part) == 1
        assert "- 1 * R_α" == part[0].to_string()

        part = r_part(order=2, exponent_of_r=0)
        assert len(part) == 1
        assert "+ 3 * R_α * R_β" == part[0].to_string()

        part = r_part(order=3, exponent_of_r=0)
        assert len(part) == 1
        assert "- 15 * R_α * R_β * R_γ" == part[0].to_string()

        part = r_part(order=4, exponent_of_r=0)
        assert len(part) == 1
        assert "+ 105 * R_α * R_β * R_γ * R_δ" == part[0].to_string()

        part = r_part(order=5, exponent_of_r=0)
        assert len(part) == 1
        assert "- 945 * R_α * R_β * R_γ * R_δ * R_ε" == part[0].to_string()

    def test_r_two_part(self):
        part = r_part(order=1, exponent_of_r=2)
        assert len(part) == 0

        part = r_part(order=2, exponent_of_r=2)
        assert len(part) == 1
        assert "- 1 * R_z^(2) * delta(α, β)" == part[0].to_string()

        part = r_part(order=3, exponent_of_r=2)
        assert len(part) == 3
        assert "+ 3 * R_z^(2) * R_α * delta(β, γ)" == part[0].to_string()

        part = r_part(order=4, exponent_of_r=2)
        assert len(part) == 6
        assert "- 15 * R_z^(2) * R_α * R_β * delta(γ, δ)" == part[0].to_string()

    def test_r_four_part(self):
        part = r_part(order=1, exponent_of_r=4)
        assert len(part) == 0
        part = r_part(order=2, exponent_of_r=4)
        assert len(part) == 0
        part = r_part(order=3, exponent_of_r=4)
        assert len(part) == 0

        part = r_part(order=4, exponent_of_r=4)
        assert len(part) == 3
        assert "+ 3 * R_z^(4) * delta(α, β) * delta(γ, δ)" == part[0].to_string()

        part = r_part(order=5, exponent_of_r=4)
        assert len(part) == 15
        assert "- 15 * R_z^(4) * R_α * delta(β, γ) * delta(δ, ε)" == part[0].to_string()


    def test_variations_between_R(self):
        result = variations_between_R_and_delta(amount_of_r=1, amount_of_delta=0)
        assert len(result) == 1
        strings = [p.to_string() for p in result]
        assert "+ 1 * R_α" in strings

        result = variations_between_R_and_delta(amount_of_r=2, amount_of_delta=0)
        assert len(result) == 1
        strings = [p.to_string() for p in result]
        assert "+ 1 * R_α * R_β" in strings

        result = variations_between_R_and_delta(amount_of_r=3, amount_of_delta=0)
        assert len(result) == 1
        strings = [p.to_string() for p in result]
        assert "+ 1 * R_α * R_β * R_γ" in strings

        result = variations_between_R_and_delta(amount_of_r=4, amount_of_delta=0)
        assert len(result) == 1
        strings = [p.to_string() for p in result]
        assert "+ 1 * R_α * R_β * R_γ * R_δ" in strings

    def test_variations_between_delta(self):
        result = variations_between_R_and_delta(amount_of_r=0, amount_of_delta=1)
        assert len(result) == 1
        strings = [p.to_string() for p in result]
        assert "+ 1 * delta(α, β)" in strings

        result = variations_between_R_and_delta(amount_of_r=0, amount_of_delta=2)
        assert len(result) == 3
        strings = [p.to_string() for p in result]
        assert "+ 1 * delta(α, β) * delta(γ, δ)" in strings
        assert "+ 1 * delta(α, γ) * delta(β, δ)" in strings
        assert "+ 1 * delta(α, δ) * delta(β, γ)" in strings

        result = variations_between_R_and_delta(amount_of_r=0, amount_of_delta=3)
        assert len(result) == 15
        strings = [p.to_string() for p in result]
        assert "+ 1 * delta(α, β) * delta(γ, δ) * delta(ε, ζ)" in strings
        assert "+ 1 * delta(α, β) * delta(γ, ζ) * delta(δ, ε)" in strings
        assert "+ 1 * delta(α, β) * delta(γ, ε) * delta(δ, ζ)" in strings

        assert "+ 1 * delta(α, γ) * delta(β, δ) * delta(ε, ζ)" in strings
        assert "+ 1 * delta(α, γ) * delta(β, ε) * delta(δ, ζ)" in strings
        assert "+ 1 * delta(α, γ) * delta(β, ζ) * delta(δ, ε)" in strings

        assert "+ 1 * delta(α, δ) * delta(β, γ) * delta(ε, ζ)" in strings
        assert "+ 1 * delta(α, δ) * delta(β, ε) * delta(γ, ζ)" in strings
        assert "+ 1 * delta(α, δ) * delta(β, ζ) * delta(γ, ε)" in strings

        assert "+ 1 * delta(α, ε) * delta(β, γ) * delta(δ, ζ)" in strings
        assert "+ 1 * delta(α, ε) * delta(β, δ) * delta(γ, ζ)" in strings
        assert "+ 1 * delta(α, ε) * delta(β, ζ) * delta(γ, δ)" in strings

        assert "+ 1 * delta(α, ζ) * delta(β, γ) * delta(δ, ε)" in strings
        assert "+ 1 * delta(α, ζ) * delta(β, δ) * delta(γ, ε)" in strings
        assert "+ 1 * delta(α, ζ) * delta(β, ε) * delta(γ, δ)" in strings


    def test_variations_between_R_and_delta(self):
        ### ONLY R ###
        self.test_variations_between_R()

        ### ONLY Delta ###
        self.test_variations_between_delta()

        ### R and Delta ###
        result = variations_between_R_and_delta(1,1)
        assert len(result) == 3
        strings = [p.to_string() for p in result]
        assert "+ 1 * R_α * delta(β, γ)" in strings
        assert "+ 1 * R_β * delta(α, γ)" in strings
        assert "+ 1 * R_γ * delta(α, β)" in strings

        result = variations_between_R_and_delta(2, 1)
        assert len(result) == 6
        strings = [p.to_string() for p in result]
        assert "+ 1 * R_α * R_β * delta(γ, δ)" in strings
        assert "+ 1 * R_γ * R_δ * delta(α, β)" in strings
        assert "+ 1 * R_β * R_γ * delta(α, δ)" in strings
        assert "+ 1 * R_β * R_δ * delta(α, γ)" in strings
        assert "+ 1 * R_α * R_γ * delta(β, δ)" in strings
        assert "+ 1 * R_α * R_δ * delta(β, γ)" in strings

        result = variations_between_R_and_delta(1, 2)
        assert len(result) == 3 * 5


    def test_construct_full_tensor(self):
        tensor = construct_all_tensor_terms(order=0)
        assert len(tensor) == 0

        tensor = construct_all_tensor_terms(order=1)
        assert len(tensor) == 1
        assert tensor[0].to_string() == "- 1 * R_α"

        tensor = construct_all_tensor_terms(order=2)
        assert len(tensor) == 2
        strings = [p.to_string() for p in tensor]
        assert "+ 3 * R_α * R_β" in strings
        assert "- 1 * R_z^(2) * delta(α, β)" in strings

        tensor = construct_all_tensor_terms(order=3)
        assert len(tensor) == 4
        strings = [p.to_string() for p in tensor]
        assert "- 15 * R_α * R_β * R_γ" in strings
        assert "+ 3 * R_z^(2) * R_α * delta(β, γ)" in strings
        assert "+ 3 * R_z^(2) * R_β * delta(α, γ)" in strings
        assert "+ 3 * R_z^(2) * R_γ * delta(α, β)" in strings

        tensor = construct_all_tensor_terms(order=4)
        assert len(tensor) == 10
        strings = [p.to_string() for p in tensor]
        assert "+ 105 * R_α * R_β * R_γ * R_δ" in strings
        assert "- 15 * R_z^(2) * R_α * R_β * delta(γ, δ)"
        assert "- 15 * R_z^(2) * R_α * R_γ * delta(β, δ)"
        assert "- 15 * R_z^(2) * R_α * R_δ * delta(β, γ)"
        assert "- 15 * R_z^(2) * R_γ * R_δ * delta(α, β)"
        assert "- 15 * R_z^(2) * R_β * R_γ * delta(α, δ)"
        assert "- 15 * R_z^(2) * R_β * R_δ * delta(α, γ)"
        assert "+ 3 * R_z^(4) * delta(α, β) * delta(γ, δ)"
        assert "+ 3 * R_z^(4) * delta(α, γ) * delta(β, δ)"
        assert "+ 3 * R_z^(4) * delta(α, δ) * delta(β, γ)"

        tensor = construct_all_tensor_terms(order=5)
        assert len(tensor) == 1+10+15

        tensor = construct_all_tensor_terms(order=6)
        assert len(tensor) == 1+15+45+15


    def test_get_product_terms_from_greek_sum(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"])], prefactor=1)
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 1
        assert result[0] == p

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "beta"])], prefactor=1)
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 1
        assert result[0] == p

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "alpha"])], prefactor=1)
        assert p.to_string() == "+ 1 * Θ_αα * Θ_αα"
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 3
        strings = [i.to_string() for i in result]
        assert "+ 1 * Θ_xx * Θ_xx" in strings
        assert "+ 1 * Θ_yy * Θ_yy" in strings
        assert "+ 1 * Θ_zz * Θ_zz" in strings

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), R("alpha"), MultipoleMoment(["alpha", "alpha"])], prefactor=2)
        p.sort_elements()
        assert p.to_string() == "+ 2 * R_α * Θ_αα * Θ_αα"
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 3
        strings = [i.to_string() for i in result]
        assert "+ 2 * R_x * Θ_xx * Θ_xx" in strings
        assert "+ 2 * R_y * Θ_yy * Θ_yy" in strings
        assert "+ 2 * R_z * Θ_zz * Θ_zz" in strings

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), R("beta"), MultipoleMoment(["alpha", "alpha"])], prefactor=1)
        p.sort_elements()
        assert p.to_string() == "+ 1 * R_β * Θ_αα * Θ_αα"
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 3
        strings = [i.to_string() for i in result]
        assert "+ 1 * R_β * Θ_xx * Θ_xx" in strings
        assert "+ 1 * R_β * Θ_yy * Θ_yy" in strings
        assert "+ 1 * R_β * Θ_zz * Θ_zz" in strings


        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "alpha"]), Delta("alpha", "beta")], prefactor=1)
        assert p.to_string() == "+ 1 * Θ_αα * Θ_αα * delta(α, β)"
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 3
        strings = [i.to_string() for i in result]
        assert "+ 1 * Θ_xx * Θ_xx * delta(x, β)" in strings
        assert "+ 1 * Θ_yy * Θ_yy * delta(y, β)" in strings
        assert "+ 1 * Θ_zz * Θ_zz * delta(z, β)" in strings


        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta"]), MultipoleMoment(["beta", "alpha"])], prefactor=1)
        assert p.to_string() == "+ 1 * Θ_αβ * Θ_αβ"
        result = get_product_terms_from_greek_sum_duplicateMultipoles(p)
        assert len(result) == 9
        strings = [i.to_string() for i in result]
        assert "+ 1 * Θ_xy * Θ_xy" in strings
        assert "+ 1 * Θ_xz * Θ_xz" in strings
        assert "+ 1 * Θ_yz * Θ_yz" in strings
        assert "+ 1 * Θ_xx * Θ_xx" in strings
        assert "+ 1 * Θ_yy * Θ_yy" in strings
        assert "+ 1 * Θ_zz * Θ_zz" in strings


    def test_get_product_terms_from_greek_sum_all(self):
        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"])], prefactor=1)
        result = get_product_terms_from_greek_sum_all(p)
        assert len(result) == 3
        strings = [i.to_string() for i in result]
        assert "+ 1 * Θ_xx" in strings
        assert "+ 1 * Θ_yy" in strings
        assert "+ 1 * Θ_zz" in strings

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha"]), MultipoleMoment(["alpha", "beta"])], prefactor=1)
        result = get_product_terms_from_greek_sum_all(p)
        assert len(result) == 3*3
        strings = [i.to_string() for i in result]
        assert '+ 1 * Θ_xx * Θ_xx' in strings
        assert '+ 1 * Θ_xx * Θ_xy' in strings
        assert '+ 1 * Θ_xx * Θ_xz' in strings#
        assert '+ 1 * Θ_xy * Θ_yy' in strings
        assert '+ 1 * Θ_yy * Θ_yy' in strings
        assert '+ 1 * Θ_yy * Θ_yz' in strings
        assert '+ 1 * Θ_xz * Θ_zz' in strings
        assert '+ 1 * Θ_yz * Θ_zz' in strings
        assert '+ 1 * Θ_zz * Θ_zz' in strings


        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "beta", "alpha"])], prefactor=1)
        assert p.to_string() == "+ 1 * Φ_ααββ"
        result = get_product_terms_from_greek_sum_all(p)
        assert len(result) == 3*3
        strings = [i.to_string() for i in result]
        assert len(set(strings)) == 3+3

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "alpha", "alpha", "alpha"])], prefactor=1)
        result = get_product_terms_from_greek_sum_all(p)
        assert len(result) == 3

        p = ProductTerm()
        p.set_elements([MultipoleMoment(["alpha", "beta", "gamma", "delta"])], prefactor=1)
        result = get_product_terms_from_greek_sum_all(p)
        assert len(result) == 3*3*3*3

