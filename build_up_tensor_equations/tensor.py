
import sympy as sp
from sympy import pretty

from Multipole import MultipoleMoment


class Tensor:
    def __init__(self, multipole:MultipoleMoment):
        self.order = len(multipole.indices)
        self.multipole = multipole
        self.indices = multipole.indices

    def get_prefactor(self):
        n = sp.symbols('n', integer=True, nonnegative=True)
        term = (-1) ** n / sp.factorial2(2 * n - 1)

        return term.subs(n, self.order)

    def get_tensor_notation(self):
        return f"T_{''.join([pretty(i.index) for i in self.indices])}"

    def to_string(self):
        return f"{self.get_prefactor()} * {self.get_tensor_notation()} * {self.multipole.to_string()}"



class TensorTerm:

    def __init__(self, m1:MultipoleMoment, m2:MultipoleMoment):
        self.tensor1 = Tensor(multipole=m1)
        self.tensor2 = Tensor(multipole=m2)

        self.order = self.tensor1.order + self.tensor2.order

        self.indices = self.tensor1.indices + self.tensor2.indices

    def get_prefactor(self):
        return self.tensor1.get_prefactor() * self.tensor1.multipole.prefactor_in_potential * self.tensor2.get_prefactor()

    def get_tensor_notation(self):
        return f"T_{''.join([pretty(i.index) for i in self.indices])}"

    def to_string(self):
        prefactor = self.get_prefactor()
        m1 = self.tensor1.multipole.to_string()
        m2 = self.tensor2.multipole.to_string()

        return f"{prefactor} * {self.get_tensor_notation()} * {m1} * {m2}"
