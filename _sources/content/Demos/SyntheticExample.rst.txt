.. _synth_example_intro:

2.1. Analysis of Simple Bodies
==============================

.. image:: ../../../Notebooks/images/SyntheticModel.png
    :align: right
    :width: 40%

As a primer to interpreting magnetic data, let’s get familiar with the
magnetic responses of some simple geologic bodies. We will then grid the
magnetic data, and investigate different visual enhancements of the data, and
apply several tools that will aid us in our geological interpretation of the
magnetic data.



.. toctree::
   :maxdepth: 1

   SyntheticMagResponse

.. image:: ./images/synthProfile.PNG
    :align: right
    :width: 40%

Here we explore the magnetic response over our synthetic geologic model. We
can view profiles through the data, and also investigate how the magnetic
response varies depending on the location of the magnetic survey and the
magnetic field orientation.


.. toctree::
   :maxdepth: 1

   SyntheticGrid

.. image:: ./images/synthGridding.png
    :align: right
    :width: 30%

This section discusses data gridding, whereby survey data are interpolated onto a regular 2D grid. This is an important initial data processing step required prior to applying filters and enhancements to magnetic data. Several gridding methods can be explored here using the data calculated from the synthetic geologic model.


.. toctree::
   :maxdepth: 1

   SyntheticVis


.. image:: ./images/synth_vis_viridis_shaded_contours.PNG
    :align: right
    :width: 30%


This section discusses image enhancement and visualization of magnetic data,
and demonstrates effects of this type of processing using the magnetic data
from the synthetic 3D geologic model. The applied filters and enhancements
reflect the initial processing steps an interpreter might take toward building
an understanding of their magnetic data and of the geology it is responding
to.


.. toctree::
   :maxdepth: 1

   SyntheticFilters

.. image:: ./images/Filters_synth_TotHoriz.PNG
    :align: right
    :width: 28%

Total field magnetic data, viewed with sun-shading or various color
enhancements, is a great approach to initially exploring a magnetic dataset.
Subtle variations within the magnetic data can be obscured however, usually by
larger or deeper magnetic bodies. Deeper or shallower sources, and more subtle
features in the magnetic data can be emphasized through the use of magnetic
data filters. This section describes and demonstrates the effect of several
commonly used magnetic data filters, including upward continuation, horizontal
and vertical derivatives, analytic signal, and tilt angle.


.. toctree::
   :maxdepth: 1

   SyntheticTiltDepth

.. image:: ./images/tilt_synthetic_edges_depths_grey.PNG
    :align: right
    :width: 33%

Interpretation of magnetic data is ideally done by geoscientists with
knowledge of the geology, lithology, and physical rock properties of typical
rock types within a project area. Commonly this is done manually, through
analysis of various magnetic data products discussed in the two previous
notebooks, and alongside other available geoscientific data. There are
however, some quick tools at our disposal to automatically pick 'edges' within
magnetic data, and help us to approximate source depth. The resulting products
can provide further guidance for geologic interpretations. This section
explains and applies one such edge and depth detection method to magnetic data
calculated from the synthetic 3D model.

