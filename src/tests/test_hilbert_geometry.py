import unittest
from hilbert import geometry
from numpy import array as na

class TestDist(unittest.TestCase):
    def test_vert(self):
        print(geometry.hdist_to_euc(na([1, 3]), na([1,5]), na([1,1]), 0.5))
    def test_horiz(self):
        print(geometry.hdist_to_euc(na([3, 1]), na([5,1]), na([1,1]), 0.5))


if __name__ == '__main__':
    unittest.main()
