.. _installation:

Installation
============

Please install the following packages to run all of the presented apps locally.
A local installation will be necessary if you wish to upload and export your own data, and if you are working with Geosoft .grd files.


Step 1: Microsoft Visual Studio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(Windows users) Update/install  `Microsoft Visual Studio 2017 <https://support.microsoft.com/en-ca/help/2977003/the-latest-supported-visual-c-downloads>`_
    - Some Python packages use C++ to speed up core functions (i.e. `SimPEG <simpeg.xyz>`_)

    .. figure:: ../images/MVS2017.PNG
        :align: center
        :width: 600


Step 2: Anaconda (Python)
^^^^^^^^^^^^^^^^^^^^^^^^^

(New Python users) `Download Anaconda <https://www.anaconda.com/download/>`_

- Launch the installation

	.. figure:: ../images/MinicondaInstaller.png
	    :align: center
	    :width: 400

- For new users, it is recommended to let Anaconda set the Environment Path

	.. figure:: ../images/AnacondaPath.png
	    :align: center
	    :width: 400


Step 3: GeoToolkit
^^^^^^^^^^^^^^^^^^

- `Download Geotoolkit <https://github.com/geoscixyz/Toolkit/archive/master.zip>`_

- Open a Command Terminal in the GeoToolkit directory (Shift+rightClick) and enter:

    .. figure:: ../images/OpenCommand.png
        :align: center
        :width: 400

- Enter>>    `conda env create -f environment.yml`

    .. figure:: ../images/InstallEnvironment.png
        :align: center
        :width: 600

Full installation time :math:`\approx 15` min

Congratulations, you should now have access to the `Python ecosystem <http://www.developintelligence.com/blog/python-ecosystem-2017/>`_!


Step 4: Run the notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^

Open a Command Terminal in the GeoToolkit directory (Shift+rightClick) and enter:

Enter>>    `jupyter notebook`

    .. figure:: ../images/LaunchNotebook.png
        :align: center
        :width: 600


You will see the list of notebooks available similar to the image below. Just click on the one you wish to work with.

Once in a notebook, you can run cells with Shift+Enter.


    .. figure:: ../images/Notebook_full_list.PNG
        :align: center
        :width: 300
        


