.. _search_example:

2.2. Case Study Using Geoscience BC SeArch II Data
==================================================

We’ve now got an idea of how some common ‘geological’ features look in magnetic data, and have been introduced to a variety of tools to help us analyze these features. Now, let’s take a look at some real magnetic data, and apply the same geophysical analysis tools that we applied to the synthetic data from our simple 3D geology scenario.

For this case study, we use a subset of data from Geoscience BC’s `SeArch Phase II`_ magnetic dataset. A study area approximately 40 km x 40 km was chosen. The magnetic data in this area may look somewhat familiar! The data from this area inspired the geologic model used in the previous case study. There is a persistent continuous magnetic high dominating the data in the east. There is a rounded magnetic feature in the center of the map area, and a linear magnetic low trending NW.


Figure of magnetic data, and geology with minfile data


Since we are dealing with ‘real’ data now, there were some ‘real life’ data preparation considerations that needed to be made before data analysis and interpretation could begin. Firstly the data needed to be cropped from the original very large dataset. We won’t actively complete data-cropping in this section, but we do provide tools to help users do this using their own data in our Toolkit gallery in Section 3. Secondly, information regarding location and dates of the magnetic survey needed to be acquired to determine the magnetic field information required for input into the processing algorithms.

Once we have the completed the processing pre-requisites, we can move forward and apply the same tools we have learned about in our previous synthetic geology case study. Again the applications are run in separate Jupyter notebooks that can be accessed through the below links:

SeArch II data case study notebooks
-----------------------------------


.. links:

.. _SeArch Phase II: http://www.geosciencebc.com/s/Report2017-03.asp