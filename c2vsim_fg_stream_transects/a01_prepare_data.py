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

# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
import os
import warnings

warnings.filterwarnings('ignore')


# Path to the root of the repository
repo_dir = os.getcwd()
current_dir = os.getcwd()

while os.path.split(repo_dir)[1] != "c2vsim_fg_stream_transects":
    repo_dir = os.path.split(repo_dir)[0]
    if os.path.split(os.path.split(repo_dir)[0])[1] == "c2vsim_fg_stream_transects":
        repo_dir = os.path.split(repo_dir)[0]

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Preparing Stream and Groundwater Data from C2VSimFG v1.5
# In this section, we prepare the model data that will be used to make the stream transects. The workflow used hereafter
# consists of the following steps:
# 1. We load the Stream Specification File and the Nodal X-Y Coordinate File and generate a lines geodataframe
# of the streams represented in the model.
# 2. Using the stream lines shapefile generated in the previous step, we select the streams of interest for our analysis.
# 3. We load the Stratigraphy File and join ground surface elevations and layer thicknesses to the stream nodes.
# 4. We load the timeseries "Groundwater Head at All Nodes" file (C2VSimFG_GW_HeadAll.out) and join simulated groundwater heads
# to stream nodes.
# 5. For each node, we extract the groundwater heads for the highest active layer. To do so, we iterate through the
# layers and select the highest one that presents a thickness greater than zero.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
