.. _installation:

Installation
============

Please install the following packages to run all of the presented apps locally.
A local installation is recommended for faster data processing, and also so that private data can be worked with locally. A local installation is also required if you plan to work with Geosoft .grd files (note: .grd file use only available to Windows users).


Step 1: Anaconda (Python)
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


Step 2: GeoToolkit
^^^^^^^^^^^^^^^^^^

- `Download Geotoolkit <https://github.com/geoscixyz/Toolkit/zipball/main/>`_

- Save the zipped file on your computer and extract the GeoToolkit folder.

- ``Double+Click`` on the ``Install_Environment.bat`` file to launch the installation:

    .. figure:: ../images/ClickInstall.png
        :align: center
        :width: 200

    .. .. figure:: ../images/InstallEnvironment.png
    ..     :align: center
    ..     :width: 600

Full installation time :math:`\approx 15` min.


Step 3: Run the notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^

- Open a Command Terminal (or 'Open PowerShell Window') in the ``Notebooks`` directory (``Shift+RightClick``) and enter:

    .. figure:: ../images/OpenCommand.png
        :align: center
        :width: 200

Enter>> ``conda activate geotoolkit``

Enter>> ``jupyter notebook``

    .. figure:: ../images/LaunchNotebook.png
        :align: center
        :width: 600


You will see the list of notebooks available similar to the image below. Just click on the one you wish to work with.

Once in a notebook, you can run cells with ``Shift+Enter``.


    .. figure:: ../images/Notebook_full_list.PNG
        :align: center
        :width: 300

Alternatively you can run the entire Notebook by selecting the ``Run All`` option from the ``Cell`` menu

    .. figure:: ../images/RunAllCells.png
        :align: center
        :width: 300
