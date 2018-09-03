import os
import testipynb
import unittest

NBDIR = os.path.sep.join(
    os.path.abspath(__file__).split(os.path.sep)[:-2] + ['docs/Notebooks']
)


class TestNotebooks(unittest.TestCase):

    def nbstripout(self):
        print('Notebook stripping out!')
        # search for notebooks
        for root, dirList, fileList in os.walk(NBDIR):
            for filename in fileList:
                if 'ipynb' in filename:
                    print('Found a ipynb to strip!')
                    os.system('nbstripout ' + os.path.join(root, filename))

    def test_notebooks(self):

        self.nbstripout()
        Test = testipynb.TestNotebooks(directory=NBDIR, timeout=10000)
        self.assertTrue(Test.run_tests())
        self.nbstripout()


if __name__ == "__main__":
    unittest.main()
