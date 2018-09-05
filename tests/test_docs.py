import subprocess
import unittest
import os

class Doc_Test(unittest.TestCase):

    @property
    def path_to_docs(self):
        dirname, filename = os.path.split(os.path.abspath(__file__))
        return os.path.sep.join(dirname.split(os.path.sep)[:-1] + ['docs'])

    def nbstripout(self):
        print('Notebook stripping out!')

        # search for images that have been missed
        for root, dirList, fileList in os.walk(self.path_to_docs):
            for filename in fileList:
                if filename.endswith(".ipynb"):
                    print('Found a ipynb to strip!')
                    os.system('nbstripout ' + os.path.join(root, filename))

    def test_html(self):
        doctrees_path = os.path.sep.join(
            self.path_to_docs.split(os.path.sep) + ['_build']+['doctrees']
        )
        html_path = os.path.sep.join(
            self.path_to_docs.split(os.path.sep) + ['_build']+['html']
        )

        self.nbstripout()

        check = subprocess.call(["sphinx-build", "-nW", "-b", "html", "-d",
            "{0!s}".format((doctrees_path)) ,
            "{0!s}".format((self.path_to_docs)),
            "{0!s}".format((html_path))])

        # self.nbstripout()

        assert check == 0



    def test_linkcheck(self):
        doctrees_path = os.path.sep.join(
            self.path_to_docs.split(os.path.sep) + ['_build']+['doctrees']
        )
        link_path = os.path.sep.join(
            self.path_to_docs.split(os.path.sep) + ['_build']
        )

        self.nbstripout()

        check = subprocess.call([
            "sphinx-build", "-nW", "-b", "linkcheck", "-d",
            "%s"%(doctrees_path),
            "%s"%(self.path_to_docs),
            "%s"%(link_path)
        ])

        # self.nbstripout()

        assert check == 0


if __name__ == '__main__':


    # print('Notebook stripping out!')
    # dirname, filename = os.path.split(os.path.abspath(__file__))
    # BASEPATH = os.path.sep.join(dirname.split(os.path.sep)[:-1])

    # # search for images that have been missed
    # for root, dirList, fileList in os.walk(self.path_to_docs):
    #     for filename in fileList:
    #         if 'ipynb' in filename:
    #             print('Found a ipynb to strip!')
    #             os.system('nbstripout ' + os.path.join(BASEPATH, filename))

    unittest.main()

    # # search for images that have been missed
    # for root, dirList, fileList in os.walk(self.path_to_docs):
    #     for filename in fileList:
    #         if 'ipynb' in filename:
    #             print('Found a ipynb to strip!')
    #             os.system('nbstripout ' + os.path.join(BASEPATH, filename))
