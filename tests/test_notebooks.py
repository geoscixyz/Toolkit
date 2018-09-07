import os
import testipynb
import unittest

NBDIR = os.path.sep.join(
    os.path.abspath(__file__).split(os.path.sep)[:-2] + ['docs/Notebooks']
)

BASEDIR = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[:-2])

def nbstripout():
    print('Notebook stripping out!')
    # search for notebooks
    for root, dirList, fileList in os.walk(BASEDIR):
        # print(root, dirList, fileList)
        for filename in fileList:
            if filename.endswith(".ipynb"):
                print('Found a ipynb to strip!')
                print(os.path.join(root, filename))
                os.system('nbstripout ' + os.path.join(root, filename))

class TestNotebooks(unittest.TestCase):

    def test_notebooks(self):
        Test = testipynb.TestNotebooks(directory=NBDIR, timeout=10000)
        self.assertTrue(Test.run_tests())
        nbstripout()


if __name__ == "__main__":
    unittest.main()
