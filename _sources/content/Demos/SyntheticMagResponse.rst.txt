.. _synth_mag_response:

2.1.1. Synthetic - Magnetic Data Response
=========================================

The geologic model
------------------

For this demonstration, a simple 3D geologic model was built attempting to capture several different types of geologic bodies.

The geologic bodies are represented by three different blocks superimposed on a background:

.. image:: ../../../Notebooks/images/SyntheticModel.png

**Block 1**: Large eastern magnetic domain, reflective of a magnetic plutonic complex or magnetic volcanic rock package.

**Block 2**: Small, 300 m\ :sup:`3` strongly magnetic domain, reflective of a shallow magnetic intrusion.

**Block 3**: North-west trending, steeply-dipping (-80 degree dip), non-magnetic feature inside the large magnetic block, reflective of a fault zone along which magnetite-destruction has occurred.

**Background**: Weakly magnetic background, reflective of weakly magnetic volcanic or sedimentary rocks.

The extents of the survey area are approximately 3 km x 3 km.




Magnetic response of a simple geologic model
--------------------------------------------

The magnetic response (total field anomaly) of the geologic model is calculated on a series of east-west flight lines running roughly perpendicular to the general strike of geologic boundaries and structures featured in the model. Flight line spacing is 200 m. Observations are at the Earth's surface. Topography is assumed to be flat in this example.

.. figure:: ./images/GeoTIFFSynthetic.png
    :align: center
    :figwidth: 40 %

An east-west profile through the magnetic data is shown below. The responses of the smaller and larger magnetic blocks are obvious in this profile, and occur directly over the sources (due to the magnetic field inclination being vertical in this example). The fault within the large eastern magnetic block is observed as a more subtle drop in the magnetic response, slightly offset from the top of the body due to its northeasterly dip.

Click below on the **'launch binder'** button or on the image to open an interactive Jupyter Notebook to further explore the magnetic response and profile of the synthetic geologic model. To launch the notebook in a new tab, right-click on the **'launch binder'** button and choose **'Open link in new tab'**. Be patient, the notebook may take a few moments to load.


.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_1_1_a_Synthetic_Mag_Data_Profile.ipynb
    :align: center

.. image:: ./images/synthEWprofile.PNG
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_1_1_a_Synthetic_Mag_Data_Profile.ipynb
    :align: center



Magnetic field effect on response
---------------------------------

As discussed in Section 1 on the Toolkit website (Magnetic Data - Background), the magnetic response will depend on the inclination, declination, and field strength of the magnetic field at the survey location. The image below shows the magnetic response at the geographic North Pole, where the Earth's magnetic field is near-vertical.

Click on the **'launch binder'** button or on the image to open an interactive Jupyter Notebook (or right-click and choose 'Open link in new tab') to investigate how the magnetic response changes with varying global locations.


.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_1_1_b_Synthetic_Mag_Data_Mag_Field.ipynb
    :align: center

.. image:: ./images/synth_mag_location_global.PNG
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_1_1_b_Synthetic_Mag_Data_Mag_Field.ipynb
    :align: center
