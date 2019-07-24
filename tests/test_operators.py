import unittest
from ..operators.operator import *

class optermTest(unittest.TestCase):
    def setUp(self):
        return

    def test_operator_addition(self):
        A = opterm(0.2, "112233")
        B = opterm(1, "112230")
        opA = operator([A])
        opB = operator([B])
        opC = opA - opB

    def test_operator_multiplication(self):
        A = opterm(1,"112233")
        B = opterm(-1, "112233")

        opA = operator([A])
        opB = operator([B])
        opC = opA * opB

    def test_scalar_multiplication(self):
        A = opterm(1,"112233")
        opA = operator([A])
        opA = 3.0*opA;
        # Assert
        opA = operator([A])
        opA = opA*3.2;
        # Assert

    def test_conjugation(self):
        A = operator([opterm(3+0.723j,"112233")])
        A = A.conj()
        self.assertEqual(A.opterms[0].string, "221133")
        self.assertEqual(A.opterms[0].coeff, 3-0.723j)

    def test_opterm_range(self):
        A = opterm(1,"000")
        self.assertEqual(A.range, 0)

        A = opterm(1,"01000")
        self.assertEqual(A.range, 1)

        A = opterm(1,"0120")
        self.assertEqual(A.range, 2)

        A = opterm(1,"112233")
        self.assertEqual(A.range, 6)

        A = opterm(1,"010003")
        self.assertEqual(A.range, 5)

if __name__ == '__main__':
  unittest.main()
