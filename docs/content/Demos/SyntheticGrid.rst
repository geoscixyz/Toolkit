.. _synth_grid:

2.1.2. Synthetic - Data Gridding
================================

Gridding the magnetic data
--------------------------

Magnetic data should be gridded prior to applying transforms and filters. This means that the flight line data must be interpolated on a regular 2D grid. There are numerous possible ways to grid data (:cite:`Briggs1974`, :cite:`Parker1983`). The minimum curvature method is shown to be robust and yield a sensible and smooth result, so we have chosen minimum curvature as the default option to complete the gridding for later filtering of the synthetic model data.

The left image below shows the data gridded using minimum curvature, and the right image shows the data gridded using nearest-neighbor.

Click on the **'launch binder'** button to launch a Jupyer notebook to try out some additional gridding methods.


.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit.git/master?filepath=.%2FNotebooks%2FSynthetic_Gridding.ipynb
    :align: center

.. image:: ./images/mincurv_nearest.png
    :target: https://mybinder.org/v2/gh/geoscixyz/Toolkit.git/master?filepath=.%2FNotebooks%2FSynthetic_Gridding.ipynb
    :align: center

