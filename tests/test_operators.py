import unittest
from ..opyrators.fermions import *

class optermTest(unittest.TestCase):
    def setUp(self):
        return

    def test_operator_addition(self):
        A = opterm(0.2, "112233")
        B = opterm(1.3, "112230")
        opA = operator([A])
        opB = operator([B])
        opC = opA + opB

        # Make sure A and B didn't change
        self.assertEqual(opA.opterms[0].string, "112233")
        self.assertEqual(opA.opterms[0].coeff, 0.2)
        self.assertEqual(opB.opterms[0].string, "112230")
        self.assertEqual(opB.opterms[0].coeff, 1.3)
        # Check resulting operator
        self.assertEqual(opC.opterms[0].coeff, 0.2)
        self.assertEqual(opC.opterms[0].string, "112233")
        self.assertEqual(opC.opterms[1].coeff, 1.3)
        self.assertEqual(opC.opterms[1].string, "112230")

    def test_identical_operator_addition(self):
        A = opterm(0.2, "112233")
        B = opterm(1.3, "112233")
        opA = operator([A])
        opB = operator([B])
        opC = opA + opB

        # Make sure A and B didn't change
        self.assertEqual(opA.opterms[0].string, "112233")
        self.assertEqual(opA.opterms[0].coeff, 0.2)
        self.assertEqual(opB.opterms[0].string, "112233")
        self.assertEqual(opB.opterms[0].coeff, 1.3)

        # Check resulting operator
        self.assertEqual(opC.opterms[0].coeff, 1.5)
        self.assertEqual(opC.opterms[0].string, "112233")

    def test_operator_subtraction(self):
        A = opterm(0.2, "112233")
        B = opterm(1.3, "112230")
        opA = operator([A])
        opB = operator([B])
        opC = opA - opB

        # Make sure A and B didn't change
        self.assertEqual(opA.opterms[0].string, "112233")
        self.assertEqual(opA.opterms[0].coeff, 0.2)
        self.assertEqual(opB.opterms[0].string, "112230")
        self.assertEqual(opB.opterms[0].coeff, 1.3)
        # Check resulting operator
        self.assertEqual(opC.opterms[0].coeff, 0.2)
        self.assertEqual(opC.opterms[0].string, "112233")
        self.assertEqual(opC.opterms[1].coeff, -1.3)
        self.assertEqual(opC.opterms[1].string, "112230")

    def test_identical_operator_subtraction(self):
        A = opterm(1.3, "112233")
        B = opterm(1.3, "112233")
        opA = operator([A])
        opB = operator([B])
        opC = opA - opB

        # Make sure A and B didn't change
        self.assertEqual(opA.opterms[0].string, "112233")
        self.assertEqual(opA.opterms[0].coeff, 1.3)
        self.assertEqual(opB.opterms[0].string, "112233")
        self.assertEqual(opB.opterms[0].coeff, 1.3)
        # Check that resulting operator is empty
        self.assertEqual(len(opC.opterms), 0)

        # Check that multiplying with a zero operator results in zero
        opD = opA * opC
        self.assertEqual(len(opD.opterms),0)

    def test_identical_operator_multiplication(self):
        A = opterm(1,"010")
        B = opterm(-1, "010")
        opA = operator([A])
        opB = operator([B])
        opC = opA * opB
        # This should result in a zero
        self.assertEqual(len(opC.opterms), 0)

    def test_cdagger_c_multiplication(self):
        A = opterm(1,"010")
        B = opterm(-1, "020")
        opA = operator([A])
        opB = operator([B])
        opC = opA * opB
        # This should result in a zero
        self.assertEqual(opC.opterms[0].string, "030")
        self.assertEqual(opC.opterms[0].coeff, -1)

    def test_c_cdagger_multiplication(self):
        A = opterm(1,"020")
        B = opterm(-1, "010")
        opA = operator([A])
        opB = operator([B])
        opC = opA * opB
        # This should result in a zero
        self.assertEqual(opC.opterms[0].string, "000")
        self.assertEqual(opC.opterms[0].coeff, -1)
        self.assertEqual(opC.opterms[1].string, "030")
        self.assertEqual(opC.opterms[1].coeff, 1)

    def test_scalar_multiplication(self):
        ##
        # Test multiplication with scalar on left
        ##
        A = opterm(1,"112233")
        opA = operator([A])
        opB = 3.0*opA;
        self.assertEqual(opA.opterms[0].string, "112233")
        self.assertEqual(opA.opterms[0].coeff, 1.0)
        self.assertEqual(opB.opterms[0].string, "112233")
        self.assertEqual(opB.opterms[0].coeff, 3.0)

        ##
        # Test multiplication with scalar on right
        ##
        opA = operator([A])
        opB = opA*3.2;
        # Assert

        self.assertEqual(opA.opterms[0].string, "112233")
        self.assertEqual(opA.opterms[0].coeff, 1.0)
        self.assertEqual(opB.opterms[0].string, "112233")
        self.assertEqual(opB.opterms[0].coeff, 3.2)

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
