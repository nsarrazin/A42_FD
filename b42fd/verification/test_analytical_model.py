import unittest
from b42fd.numerical_model.model import *
from b42fd.verification.analytical_model import *

class TestAnalytical(unittest.TestCase):


    def test_calculate_eigenvalues_1(self):
        """hand-verifiable test cases are used here"""
        A = 1
        B = 2
        C = 1

        V0 = 2
        c = 1
        hand_real = -2
        hand_imag = 0
        res = calculate_eigenvalues(A, B, C, V0, c)
        self.assertEqual(res[0], complex(hand_real, hand_imag))
        self.assertEqual(res[1], complex(hand_real, -1*hand_imag))

    def test_calculate_eigenvalues_2(self):
        """hand-verifiable test cases are used here"""
        A = 1
        B = 4
        C = 5

        V0 = 3
        c = 2
        hand_real = -3
        hand_imag = 1.5
        res = calculate_eigenvalues(A, B, C, V0, c)
        self.assertEqual(res[0], complex(hand_real, hand_imag))
        self.assertEqual(res[1], complex(hand_real, -1*hand_imag))


    def test_short_period_eigenvalues_1(self):

        """these values are hand-verified. the calculate_eigenvalues function is assumed
        to be verified in a previous unit test and so it's used here. integers are desired for the hand results so the
        result can be compared exactly. the assertion breaks if the values are not in the same order, which is another
        reason the calculate_eigenvalues function is preferred"""
        test_data = {'muc': 1000, 'KY2': 2, 'CZadot': -.001, 'CZa': -5, 'CZq': -5, 'Cmq': -10, 'Cma': -1, 'V0': 1,
                     'c': 1}
        hand_A = 8000004
        hand_B = 18005
        hand_C = 1895
        hand_result = calculate_eigenvalues(hand_A, hand_B, hand_C, test_data['V0'], test_data['c'])
        self.assertEqual(short_period_eigenvalues(test_data), hand_result)

    def test_short_period_eigenvalues_2(self):
        """these values are hand-verified. the calculate_eigenvalues function is assumed
        to be verified in a previous unit test and so it's used here. integers are desired for the hand results so the
        result can be compared exactly. the assertion breaks if the values are not in the same order, which is another
        reason the calculate_eigenvalues function is preferred"""
        test_data = {'muc': 1500, 'KY2': 3, 'CZadot': -.002, 'CZa': -4, 'CZq': -4, 'Cmq': -10, 'Cma': -1, 'V0': 4,
                     'c': 1}
        hand_A = 27000018
        hand_B = 33004
        hand_C = 2916
        hand_result = calculate_eigenvalues(hand_A, hand_B, hand_C, test_data['V0'], test_data['c'])
        self.assertEqual(short_period_eigenvalues(test_data), hand_result)

    def test_phugoid_1(self):
        test_data = {'muc': 100, 'CZa': -5, 'Cmq': -10, 'Cma': -1, 'CXu': -1, 'CZu': -1, 'CXa': 1, 'CZ0': -10, 'Cmu': 1,
                     'V0': 2, 'c': 1}
        hand_A = 50000
        hand_B = 60
        hand_C = 60
        hand_result = calculate_eigenvalues(hand_A, hand_B, hand_C, test_data['V0'], test_data['c'])

        self.assertEqual(phugoid_eigenvalues(test_data), hand_result)

    def test_phugoid_2(self):
        test_data = {'muc': 10, 'CZa': -10, 'Cmq': -10, 'Cma': -2, 'CXu': -1, 'CZu': -3, 'CXa': 1, 'CZ0': -100,
                     'Cmu': 10, 'V0': 10, 'c': 1}
        hand_A = 2800
        hand_B = -30
        hand_C = 10600
        hand_result = calculate_eigenvalues(hand_A, hand_B, hand_C, test_data['V0'], test_data['c'])

        self.assertEqual(phugoid_eigenvalues(test_data), hand_result)

    def test_dutchroll_1(self):
        test_data = {'mub': 100, 'KX2': 2, 'Cnr': -1, 'Cnb': 1, 'V0': 10, 'c': 2}
        hand_A = -400
        hand_B = -.5
        hand_C = -1
        hand_result = calculate_eigenvalues(hand_A, hand_B, hand_C, test_data['V0'], test_data['c'])

        self.assertEqual(dutch_roll_eigenvalues(test_data), hand_result)

    def test_dutchroll_2(self):
        test_data = {'mub': 400, 'KX2': 7, 'Cnr': -1, 'Cnb': 2, 'V0': 100, 'c': 5}
        hand_A = -5600
        hand_B = -.5
        hand_C = -2
        hand_result = calculate_eigenvalues(hand_A, hand_B, hand_C, test_data['V0'], test_data['c'])

        self.assertEqual(dutch_roll_eigenvalues(test_data), hand_result)



if __name__ == '__main__':
    unittest.main()
