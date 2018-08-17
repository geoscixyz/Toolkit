import sys
import os

def nbstripout():
    # get relevant directories
    cwd = os.getcwd()

    # search for images that have been missed
    for root, dirList, fileList in os.walk(cwd):
        for filename in fileList:
            if 'ipynb' in filename:
                os.system('nbstripout ' + os.path.join(root, filename))


nbstripout()
