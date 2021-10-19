.. _search_grid_vis:

2.2.1. Search Case Study - Magnetic Data Visualization
======================================================


Magnetic response over the Search Phase II project area
-------------------------------------------------------

Gridded residual magnetic data from Geoscience BC's Search Phase II project (Anomalous Magnetic Field) were downloaded from the Search II project page on Geoscience BC's website. The magnetic data is interpolated onto a 50 m\ :sup:`2` grid. The data are shown using coordinate system NAD 83 UTM Zone 9.

.. figure:: ./images/SearchQuestII.png
    :align: center
    :figwidth: 100 %


Click below on the **'launch binder'** button or on the image to access an interactive notebook allowing you to explore application of various visual enhancements to the Search Phase II data. Continue below for some examples of image enhancement as applied to the Search II data.


.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_2_1_Search_Mag_Data_Visualization.ipynb
    :align: center

.. image:: ./images/search_vis_notebook_snapshot.PNG
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_2_1_Search_Mag_Data_Visualization.ipynb
    :align: center


The interactive notebook allows users to choose a subset of magnetic data from the larger Search II dataset.


.. figure:: ./images/Search_subset_and_TMI.PNG
    :align: center
    :figwidth: 100 %


.. figure:: ./images/search_subset_profile.PNG
    :align: center
    :figwidth: 100 %


Data visualization and image enhancement
----------------------------------------

As in the synthetic model example presented previously (:ref:`Section 2.1.3<synth_filters>`), we can explore different ways of presenting and enhancing the magnetic image to find optimal visual parameters for interpreting features of interest.

For example, here is the specified subset of data with and without sun-shading applied:

.. figure:: ./images/search_vis_without_and_with_shading.PNG
    :align: center
    :figwidth: 100 %


Here is the data with red-blue (upper left), viridis (upper right), jet (lower left), and grey scale (lower right) colour maps applied

.. figure:: ./images/search_vis_colour_maps.png
    :align: center
    :figwidth: 100 %

Again, either a linear or a histogram equalized colour stretch can be applied:

.. figure:: ./images/search_stretch_linear_histoeq.png
    :align: center
    :figwidth: 100 %

**The interactive notebook will also allow you to export maps with your chosen colour maps and visual enhancements**. These maps are exported as located GeoTiff files and can be loaded into any software that accepts these types of files, including Google Earth Pro and ArcMap. Your exported maps will be stored on the cloud temporarily if you are running the notebooks online. You will need to install the Toolkit locally to enable local saving of exported maps. Please follow this link and Toolkit installation instructions :ref:`here<installation>`

The below image shows a located image of the Search II magnetic data in Google Earth Pro.

.. figure:: ./images/draped_mag_search.PNG
    :align: center
    :figwidth: 100 %
