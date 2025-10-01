# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: editable,slideshow,tags,-all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
# ---

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Discussion
# Our extracted streambed elevations from bathymetry show significant spatial variability, likely because stream nodes
# do not perfectly align with the lowest points of the channel. A potential improvement would be refining the transect
# using bathymetry-informed adjustments or applying a DEM-based stream delineation algorithm to connect the lowest points
# of the channel along the flow direction.
#
# Our approach to classifying gaining vs. losing reaches by comparing groundwater levels to bathymetry is reasonable for
# shallow sections where bathymetry approximates river stage. However, in deeper, low-slope sections downstream, the
# difference between streambed elevation and actual water stage can be substantial. Incorporating interpolated stage
# time series from stream gages in future iterations would provide an additional benchmark.
#
