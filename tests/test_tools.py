import unittest
from misc import tools

class MyTestCase(unittest.TestCase):
    def test_binary_search_rotate(self):

        s = tools.BinarySearcher(7, 19, discrete=True, rotate=25)
        # find '9'
        n = -1
        while n != 9:
            n = s.next()
            s.feedback(higher=n<9)
            print(n)

        print("found! ", n)







































if __name__ == '__main__':
    unittest.main()
