import subtiwiki
import unittest

class SubtiwikiTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_hypergeometric(self):
        results = subtiwiki.hypergeometric_test(['xkdK', 'xkdG', 'xkdF', 'xkdO'], return_header=False)
        print(results)
        self.assertTrue(len(results) > 0)
        self.assertTrue(results[0][0] == 'Xpf Regulon sigma factor')

if __name__ == '__main__':
    unittest.main()
