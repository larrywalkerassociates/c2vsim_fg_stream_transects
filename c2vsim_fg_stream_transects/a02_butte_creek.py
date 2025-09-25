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

import datetime
import json
import os
import requests
import warnings

from shapely import LineString, Point
import geopandas as gpd
import pandas as pd

from iwfm_lwa import (skip_until_flag, get_groundwater_nodes, load_gwalloutfl, get_stratigraphy, get_stream_nodes,
                      get_var, line_to_list)

warnings.filterwarnings('ignore')


# Path to the root of the repository
repo_dir = os.getcwd()
current_dir = os.getcwd()

while os.path.split(repo_dir)[1] != "c2vsim_fg_stream_transects":
    repo_dir = os.path.split(repo_dir)[0]
    if os.path.split(os.path.split(repo_dir)[0])[1] == "c2vsim_fg_stream_transects":
        repo_dir = os.path.split(repo_dir)[0]

os.chdir(repo_dir)

from functions.make_stream_lines_gdf import make_stream_lines_gdf
from functions.project_points_gdf_to_line_string import project_points_gdf_to_line_string

os.chdir(current_dir)

# Let's set key paths
data_dir = os.path.join(repo_dir, "data")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Butte Creek Stream Transect
# In this section, we will make a stream transect for Butte Creek using the C2VSim-FG model streambed and groundwater
# table elevations, . The transect will be used to
# In this section, we prepare the model data that will be used to make the stream transects. The workflow used hereafter
# consists of the following steps:
# 1. We load the Stream Specification File and the Nodal X-Y Coordinate File and generate a lines geodataframe
# of the streams represented in the model.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
