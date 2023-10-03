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




    def test_binary_search(self):
        for end in range(2, 20):
            for start in range(0, 10):
                for treasure in range(start, end):
                    s = tools.BinarySearcher(start, end, discrete=True)
                    while s.has_next():
                        n = s.next()
                        if n == treasure:
                            break
                        s.feedback(higher=n<treasure)





































if __name__ == '__main__':
    unittest.main()
