.. _synth_filters:

2.1.3. Synthetic - 2D Magnetic Data Filters
===========================================

2D magnetic data filters are commonly applied to magnetic data to allow us to view that data in different ways. They can highlight shallower or deeper magnetic sources, and can emphasize gradients in the magnetic data that occur at geological contacts and in association with structure.

We use the same synthetic model here as in the previous notebook to demonstrate the effects of various 2D magnetic filters.

Short descriptions of select filters are provided in the sections below. Summaries of various filters and their role in gravity and magnetic data interpretation can be found in Roest et al., 1992, :cite:`Macleod1993`, Miller and Singh, 1994, Milligan and Gunn, 1997, Verduzco et al., 2004, Isles and Rankin, 2013, Dentith and Mudge, 2014.  

Click here to go to the interactive magnetic data filtering notebook, and try applying the various filters discussed below.

Upward continuation
-------------------

Upward continuation can be considered a ‘wavelength’ filter. Upward continuation simulates the magnetic response that would be observed if data were collected at a greater height above the Earth’s surface than it was actually collected. The result is that longer-wavelength, deeper features are emphasized over shorter-wavelength, near-surface features. It is effective for interpretation of deeper geology, or for reducing noise that may be found in data collected very close to the ground. 



X and Y derivatives of the magnetic response
--------------------------------------------

Calculating horizontal and vertical derivatives of the magnetic response gives us a way to better visualize local gradients in the data that can be obscured by larger and deeper features in the total field response. X and Y derivatives are calculated based on the difference observed in the magnetic response between adjacent points, or grid cells, in the X and Y directions, respectively: 


.. math::

	\frac{\partial B}{\partial x}(x_n,y_n)= \frac{B_{(x_{n+1},y_n)}-B_{(x_{n-1},y_n)}}{2\Delta x}


.. math::

	\frac{\partial B}{\partial y}(x_n,y_n)= \frac{B_{(x_n,y_{n+1})}-B_{(x_n,y_{n-1})}}{2\Delta y}


where B is the magnetic field. Steep gradients usually occur at geologic contacts, faults, or fractures, where there is a distinct contrast in magnetic susceptibility. Interpreters should not consider X and Y derivative results independently, and should be sure to evaluate both, since visualizing gradients in one direction will de-emphasize features trending in the perpendicular direction. The X and Y derivatives are also used in the calculation of several other magnetic data filters.


Vertical derivative
-------------------

A very commonly used magnetic filter is the first vertical derivative. This process considers the change in the magnetic response when calculated at two different heights above the ground. 

The effect is to emphasize near-surface features where the difference between the responses calculated at two different heights will be largest. First vertical derivative anomalies will highlight edges of magnetic sources, and appear over the top of the source if the contact or feature is vertical. 


Total horizontal derivative
---------------------------

The total horizontal derivative combines the X and Y derivatives: 

.. math::

	\frac{\partial B}{\partial r}(x_n,y_n) = \sqrt{\left({\frac{\partial B}{\partial x}}\right)^2 + \left({\frac{\partial B}{\partial y}}\right)^2}

Gradients in both directions are now accounted for in a single map. The highest total horizontal derivative values occur at the edges or boundaries of magnetic sources. It does not detect narrow sources as effectively as the vertical derivative.

Tilt angle
----------


The tilt derivative, or tilt angle, normalizes the vertical derivative by the horizontal derivatives: 


.. math::

	TDR(x,y) = tan^{-1}\left[\frac{\partial B}{\partial z}\Bigg/{\sqrt{\left({\frac{\partial B}{\partial x}}\right)^2 + \left({\frac{\partial B}{\partial y}}\right)^2}}\right]


A magnetic source exhibiting a strong contrast with surrounding rocks (e.g. a large, near surface, magnetic unit) will yield both high vertical and horizontal gradients, and a more weakly contrasting body will yield proportionally smaller vertical and horizontal gradients. Normalizing the vertical by horizontal derivatives means different amplitude responses are assigned equivalent values. This is a very useful filter for enhancing more subtle features in the magnetic data. The tilt angle ranges from -90 to +90 degrees, is positive over the magnetic source, and negative outside the source, with the body’s edge delimited by the 0 degree contour.

Analytic signal
---------------

Analytic signal is also known as the total gradient. It is calculated from the vertical and horizontal derivatives: 


.. math::

	AS(x,y) = \sqrt{\left({\frac{\partial B}{\partial x}}\right)^2 + \left({\frac{\partial B}{\partial y}}\right)^2 + \left({\frac{\partial B}{\partial z}}\right)^2}


The advantage is that we now capture gradients in all three directions, further enhancing detection of geologic boundaries and structures. The analytic signal peaks above narrow bodies and along the edges of larger geologic features that are in magnetic contrast to their surroundings. The analytic signal is a very useful interpretation product in areas where magnetic remanence is suspected, and in areas of low latitude since it is not affected by magnetisation direction.