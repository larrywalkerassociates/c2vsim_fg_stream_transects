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
# # Introduction
#  Over the past decade, SGMA implementation in California has generated vast amounts of observed and simulated data.
#  Presenting these data in a meaningful, intuitive, and transparent way requires advanced processing and visualization
#  tools. However, most mainstream tools are proprietary, leading to recurring license costs and limited flexibility.
#
#  The widespread adoption of high-level programming such as Python and R has enabled the development of open-source
#  libraries that support customizable workflows—from data processing and visualization to interactive dashboards.
#
#  In this work, we present a Python-based workflow to integrate public datasets of observed and simulated groundwater
#  levels, lithology logs, and bathymetry, to create animated stream transects for the Sacramento River. These animations
#  help assess stream-groundwater connectivity under varied hydrologic conditions varied lithologies. The workflow is
#  customizable and replicable for other streams across California's Central Valley.
