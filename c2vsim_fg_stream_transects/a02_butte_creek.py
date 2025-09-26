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

import functools
import os
import warnings

import geopandas as gpd
import ipywidgets as widgets
import pandas as pd


warnings.filterwarnings('ignore')


# Path to the root of the repository
repo_dir = os.getcwd()
current_dir = os.getcwd()



while os.path.split(os.path.split(repo_dir)[0])[1] == "c2vsim_fg_stream_transects":
    repo_dir = os.path.split(repo_dir)[0]


os.chdir(repo_dir)

from functions.make_transect_animation import make_transect_animation
from functions.make_transect_layout import make_transect_layout

os.chdir(current_dir)

# Let's set key paths
data_dir = os.path.join(repo_dir, "data")


# Path to Butte Creek stream nodes file
butte_creek_stream_nodes_csv_path = os.path.join(data_dir, "butte_creek_stream_nodes.csv")

# Path to Butte Creek simulated groundwater heads
butte_creek_wt_stream_nodes_csv_path = os.path.join(data_dir, "butte_creek_gwallout_stream_nodes.csv")

# Path to Butte Creek CASGEM observations lookup table
obs_wells_butte_creek_shp_path = os.path.join(data_dir, "obs_wells_butte_creek.shp")

# Path to Butte Creek CASGEM observed groundwater levels
obs_wells_butte_creek_csv_path = os.path.join(data_dir, "obs_wells_butte_creek.csv")

# Path to Butte Creek lithology logs shapefile
lith_logs_butte_creek_shp_path = os.path.join(data_dir, "lith_logs_butte_creek.shp")

# Path to Butte Creek AEM lithology logs
butte_lith_csv_path = os.path.join(data_dir, "butte_creek_lithology.csv")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Butte Creek Stream Transect
# In this section, we will make a stream transect for Butte Creek using the C2VSim-FG model streambed and groundwater
# table elevations, DWR CASGEM groundwater levels observations, and DWR's AEM lithology logs.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]

# Let's load the stream nodes for Butte Creek
stream_nodes_df = pd.read_csv(butte_creek_stream_nodes_csv_path)

# Let's load the simulated groundwater table at the stream nodes
sim_df = pd.read_csv(
    butte_creek_wt_stream_nodes_csv_path,
    parse_dates=["date"]
)

# Let's load the lithology logs
lith_df = pd.read_csv(butte_lith_csv_path)

# Let's load the lithology logs shapefile
lith_lut_df = gpd.read_file(lith_logs_butte_creek_shp_path)

obs_lut_df = gpd.read_file(obs_wells_butte_creek_shp_path)

obs_df = pd.read_csv(obs_wells_butte_creek_csv_path, parse_dates=["date"])


kwargs = {}
html_video = make_transect_animation(
        stream_nodes_df,
        sim_df,
        obs_df,
        obs_lut_df,
        lith_lut_df,
        lith_df,
    **kwargs)

widgets.HTML(
        value=html_video,
        placeholder="",
        description="",
    )
