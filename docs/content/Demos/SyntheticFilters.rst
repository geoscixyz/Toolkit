.. _synth_filters:

2.1.4. Synthetic - 2D Magnetic Data Filters
===========================================

2D magnetic data filters are commonly applied to magnetic data to allow us to view that data in different ways. They can highlight shallower or deeper magnetic sources, and can emphasize gradients in the magnetic data that occur at geological contacts and in association with structure.

We use the same synthetic model here as in the previous notebook to demonstrate the effects of various 2D magnetic filters.

Short descriptions of select filters are provided in the sections below. Summaries of various filters and their role in gravity and magnetic data interpretation can be found in :cite:`Nabighian1972`, :cite:`Roest1992`, :cite:`Macleod1993`, :cite:`Miller1994`, :cite:`Milligan1997`, :cite:`Verduzco2004`, :cite:`Salem2005`, :cite:`Isles2013`, :cite:`Dentith2014`.

Click on the **'launch binder'** button below (or right-click and choose 'Open link in new tab') to go to the interactive magnetic data filtering notebook, and try applying the various filters discussed below.

.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_1_4_Synthetic_Mag_Data_Filters.ipynb
    :align: center

.. image:: ./images/synth_filter_notebook_snapshot.PNG
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit/main?filepath=.%2FNotebooks%2F2_1_4_Synthetic_Mag_Data_Filters.ipynb
    :align: center


Upward continuation
-------------------

Upward continuation can be considered a ‘wavelength’ filter. **Upward continuation simulates the magnetic response that would be observed if data were collected at a greater height above the Earth’s surface than it was actually collected**. The result is that longer-wavelength, deeper features are emphasized over shorter-wavelength, near-surface features. It is effective for interpretation of deeper geology, or for reducing noise that may be found in data collected very close to the ground.

In our synthetic example, the small magnetic anomaly and the narrow fault feature quickly become indistinct with increased survey height (upward continued to 50 m), however the larger deeper magnetic body persists in the data. The left image below shows the TMI, and the right image shows the upward continued magnetic data.

.. figure:: ./images/Filters_synth_TMI_upcontin50.png
    :align: center
    :figwidth: 100 %


.. _synthfilters_XY_deriv:

X and Y derivatives of the magnetic response
--------------------------------------------

Calculating horizontal and vertical derivatives of the magnetic response gives us a way to better visualize local gradients in the data that can be obscured by larger and deeper features in the total field response. **X and Y derivatives are calculated based on the difference observed in the magnetic response between adjacent points**, or grid cells, in the X and Y directions, respectively:


.. math::

	\frac{\partial B}{\partial x}(x_n,y_n)= \frac{B_{(x_{n+1},y_n)}-B_{(x_{n-1},y_n)}}{2\Delta x}


.. math::

	\frac{\partial B}{\partial y}(x_n,y_n)= \frac{B_{(x_n,y_{n+1})}-B_{(x_n,y_{n-1})}}{2\Delta y}


where B is the magnetic field. **Steep gradients usually occur at geologic contacts, faults, or fractures, where there is a distinct contrast in magnetic susceptibility**. Interpreters should not consider X and Y derivative results independently, and should be sure to evaluate both, since visualizing gradients in one direction will de-emphasize features trending in the perpendicular direction. The X and Y derivatives are also used in the calculation of several other magnetic data filters.

Grids showing the X (left image) and Y (right image) derivatives are shown below. Note the emphases on east and west geological boundaries using the X derivative, and on the north and south geological boundaries using the Y derivative (best seen in relation to the smaller magnetic anomaly).


.. figure:: ./images/Filters_synth_XY_Deriv.png
    :align: center
    :figwidth: 100 %


.. _synthfilters_vert_deriv:

Vertical derivative
-------------------

A very commonly used magnetic filter is the first vertical derivative. This process **considers the change in the magnetic response when calculated at two different heights above the ground**.

The effect is to emphasize near-surface features where the difference between the responses calculated at two different heights will be largest. **First vertical derivative anomalies will highlight edges of magnetic sources, and appear over the top of the source if the contact or feature is vertical**.

Notice how the first vertical derivative peaks directly over the small magnetic body from the synthetic geologic model. It also highlights the edge, and main body, of the large eastern magnetic block. The de-magnetized NW-trending fault inside the eastern block appears offset in the data as a result of it's north-easterly dip.

.. figure:: ./images/Filters_synth_VertDeriv.PNG
    :align: center
    :figwidth: 50 %

.. _synthfilters_tot_horiz_deriv:

Total horizontal derivative
---------------------------

The total horizontal derivative combines the X and Y derivatives:

.. math::

	\frac{\partial B}{\partial r}(x_n,y_n) = \sqrt{\left({\frac{\partial B}{\partial x}}\right)^2 + \left({\frac{\partial B}{\partial y}}\right)^2}

Gradients in both directions are now accounted for in a single map. **The highest total horizontal derivative values occur at the edges or boundaries of magnetic sources**. It does not detect narrow sources as effectively as the vertical derivative.

The total horizonal derivative can be seen to peak here over the edges of the small magnetic body, and along the margin of the large magnetic block in the east. Again, the magnetic patterns internal to the eastern block are complicated due to the dipping nature of the demagnetized fault.

.. figure:: ./images/Filters_synth_TotHoriz.PNG
    :align: center
    :figwidth: 50 %


.. _synthfilters_tilt_angle:

Tilt angle
----------


The tilt derivative, or tilt angle, **normalizes the vertical derivative by the horizontal derivatives**:


.. math::

	TDR(x,y) = tan^{-1}\left[\frac{\partial B}{\partial z}\Bigg/{\sqrt{\left({\frac{\partial B}{\partial x}}\right)^2 + \left({\frac{\partial B}{\partial y}}\right)^2}}\right]


A magnetic source exhibiting a strong contrast with surrounding rocks (e.g. a large, near surface, magnetic unit) will yield both high vertical and horizontal gradients, and a more weakly contrasting body will yield proportionally smaller vertical and horizontal gradients. **Normalizing the vertical by horizontal derivatives means different amplitude responses are assigned equivalent values. This is a very useful filter for enhancing more subtle features in the magnetic data**. The tilt angle ranges from -90 to +90 degrees, is **positive over the magnetic source, and negative outside the source, with the body’s edge delimited by the 0 degree contour**.

It is helpful to view the tilt angle using a color map that highlights the middle, or near-zero tilt angle values, which should trace source edges. For example, the red-blue color map (RdBu).

In the grid image below, you can see that zero values trace the edge of the small magnetic body, with positive tilt angle values directly over the top of (inside) the body and negative values outside the magnetic body. Zero values also mark the edge of the large eastern magnetic block.

.. figure:: ./images/Filters_synth_TiltAngle.PNG
    :align: center
    :figwidth: 50 %


.. _synthfilters_an_sig:

Analytic signal
---------------

Analytic signal is also known as the total gradient. It is calculated from the vertical and horizontal derivatives:


.. math::

	AS(x,y) = \sqrt{\left({\frac{\partial B}{\partial x}}\right)^2 + \left({\frac{\partial B}{\partial y}}\right)^2 + \left({\frac{\partial B}{\partial z}}\right)^2}


The advantage is that we now capture gradients in all three directions, further enhancing detection of geologic boundaries and structures. **The analytic signal peaks above narrow bodies and along the edges of larger geologic features that are in magnetic contrast to their surroundings**. **The analytic signal is a very useful interpretation product in areas where magnetic remanence is suspected, and in areas of low latitude since it is not affected by magnetisation direction**.

The analytic signal from the synthetic model looks similar to the first vertical derivative since an anomaly is mapped directly over the smaller magnetic body, however it differs in that it maps the edge of the larger eastern magnetic block but not the top.


.. figure:: ./images/Filters_synth_AnSig.PNG
    :align: center
    :figwidth: 50 %
