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
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from myst_nb import glue
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
from functions.add_bore_logs import add_bore_logs

os.chdir(current_dir)

# Let's set key paths
data_dir = os.path.join(repo_dir, "data")

stream_name = "SACRAMENTO"

glue("stream_name", stream_name.capitalize(), display=False)

# Path to stream nodes file
stream_nodes_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_stream_nodes.csv")

# Path to simulated groundwater heads
wt_stream_nodes_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_gwallout_stream_nodes.csv")

# Path to CASGEM observations lookup table
obs_wells_shp_path = os.path.join(data_dir, f"obs_wells_{stream_name.lower()}.shp")

# Path to CASGEM observed groundwater levels
obs_wells_csv_path = os.path.join(data_dir, f"obs_wells_{stream_name.lower()}.csv")

# Path to lithology logs shapefile
lith_logs_shp_path = os.path.join(data_dir, f"lith_logs_{stream_name.lower()}.shp")

# Path to AEM lithology logs
lith_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_lithology.csv")


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Southern Transect (Sacramento to Yuba City)
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]

# Let's load the stream nodes
stream_nodes_df = pd.read_csv(stream_nodes_csv_path)

# Let's load the simulated groundwater table at the stream nodes
sim_df = pd.read_csv(
    wt_stream_nodes_csv_path,
    parse_dates=["date"]
)

# Let's load the lithology logs
lith_df = pd.read_csv(lith_csv_path)

# Let's load the lithology logs shapefile
lith_lut_df = gpd.read_file(lith_logs_shp_path)

obs_lut_df = gpd.read_file(obs_wells_shp_path)

# Let's remove deep piezometers of nested well 386062N1215602 and adjacent deep wells
# 386121N1215550W001, 386008N1215488W001, and 386061N1215313W001 and in crowded area
# north of Sacramento

obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "386062N1215602W002",
        "386121N1215550W001",
        "386061N1215313W001",
        "386008N1215488W001"
    ]
)].reset_index(drop=True)

# Let's remove deep wells 386449N1215631W002, 386489N1215679W001, and 386494N1215650W001
# near Natomas

obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "386449N1215631W002",
        "386489N1215679W001",
        "386494N1215650W001"
    ]
)].reset_index(drop=True)

# Let's remove well without depth information south of I5
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "386600N1216157W001"
    ]
)].reset_index(drop=True)

# Let's remove deep wells north of I-5
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "386774N1216353W001",
        "386790N1216348W001",
        "386821N1216345W002"
    ]
)].reset_index(drop=True)
# Let's move a bit the projected point of one of the monitoring
# wells north of I5 to prevent overlapping of labels
obs_lut_df.loc[
    obs_lut_df["site_code"]=="386772N1216337W001",
    "proj"
] = 355000

obs_lut_df.loc[
    obs_lut_df["site_code"]=="386694N1216153W002",
    "proj"
] = 357000

# Let's remove irrigation well by SMF
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "387033N1216181W002"
    ]
)].reset_index(drop=True)

# Let's remove deepest wells by Fremont weir
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "387552N1215951W001",
        "387566N1215947W001",
        "387605N1215946W001"
    ]
)].reset_index(drop=True)

# Remove irrigation wells by Sutter bypass
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "387725N1216004W001",
        "387783N1216078W001"
    ]
)].reset_index(drop=True)

# Let's remove too shallow piezometer for nested wells 387626N1216357,
# 387628N1216350
# and irrigation wells by Fremont Weir

obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "387626N1216357W001",
        "387658N1216311W001",
        "387721N1216311W001",
        "387628N1216350W001",
        "387630N1216336W001"
    ]
)].reset_index(drop=True)

obs_lut_df.loc[
    obs_lut_df["site_code"]=="387626N1216357W002",
    "proj"
] = 338500

obs_lut_df.loc[
    obs_lut_df["site_code"]=="387630N1216336W002",
    "proj"
] = 337000

obs_lut_df.loc[
    obs_lut_df["site_code"]=="387628N1216350W002",
    "proj"
] = 335500

# Let's remove a deep well by Kings Landing
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "388010N1217279W001"
    ]
)].reset_index(drop=True)

# Let's remove some deep wells north of Kings Landing
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "388373N1217344W001",
        "388375N1217343W001",
        "388393N1217330W001",
        "388393N1217330W003",
        "388393N1217330W004"
    ]
)].reset_index(drop=True)

# Let's remove deep piezometers of nested well by State Ranch Bend
# and adjacent deep well
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "388813N1217525W003",
        "388813N1217525W002",
        "388730N1217474W001"
    ]
)].reset_index(drop=True)

# Let's remove well without depth data south of Tyndall Landing
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "388643N1218042W001"
    ]
)].reset_index(drop=True)

# Let's remove well without depth data north of Kirkville
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "389191N1218102W001"
    ]
)].reset_index(drop=True)

# Let's remove deepest well by Poffenbergers Landing
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "389382N1218291W001"
    ]
)].reset_index(drop=True)

# Let's remove deepest well by Millers Landing
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "389596N1218314W001"
    ]
)].reset_index(drop=True)

# Let's remove deep wells by Fraziers landing
obs_lut_df = obs_lut_df.loc[~obs_lut_df["site_code"].isin(
    [
        "390113N1218323W001",
        "390121N1218304W001",
        "390124N1218291W001",
        "390124N1218291W002",
        "390124N1218291W003"
    ]
)].reset_index(drop=True)

# Let's move a well that is slightly out of the margins
obs_lut_df.loc[
    obs_lut_df["site_code"]=="390587N1218380W001",
    "proj"
] = 277000

# Let's load the observed groundwater levels

obs_df = pd.read_csv(obs_wells_csv_path, parse_dates=["date"])

obs_df = obs_df.loc[obs_df["site_code"].isin(obs_lut_df["site_code"])].reset_index(drop=True)


step_miles = 1 * np.power(
    10, np.modf
            (np.log10(
            (
                    stream_nodes_df["proj"].max() - stream_nodes_df["proj"].min()
            ) / 1609.34
        )
        )[1] - 1
)

transect_start = stream_nodes_df["proj"].max()-50*1609.34  # 120 miles from south end, Sacramento

transect_end = stream_nodes_df["proj"].max()-120*1609.34 # 190 miles from south end, Yuba

transect_length = transect_start - transect_end

ticks = np.arange(
    0,
    transect_length,
    step_miles * 1609.34)
labels = [f"{int((stream_nodes_df['proj'].max()-transect_start + tick) / 1609.34):,.0f}" for tick in ticks]
xticks_kwargs = {
    "ticks": transect_start - ticks,
    "labels": labels
}

tick_params_kwargs = {
    "labelsize": "xx-large"
}

xlabel_kwargs = {
    "fontsize": "xx-large"
}

ylabel_kwargs = {
    "fontsize": "xx-large"
}

title_kwargs = {
    "pad": 150,
    "fontsize": "xx-large"
}

legend_kwargs = {
                "loc": "lower center",
                "bbox_to_anchor": (0.5, 0),
                "ncols": 4,
    "fontsize": "xx-large"
            }

sim_kwargs = {
    "color": "blue",
    "label": "Simulated GWE",
    "linewidth": 2
}

subplot_kwargs = {"nrows": 2,
                              "height_ratios": [0.7, 0.3],
                              "figsize": (12, 8)}

obs_kwargs = {
            "label": 'Observed GWE',
            "markerfacecolor": "#4B7164",
            "marker": 'o',
            "color": "w",
    "markersize": 5
        }

streambed_kwargs = {
            "color": "#555555",
            "linestyle": "dashed",
            "label": "Bathymetry"
        }

kwargs = {
    "tick_params_kwargs": tick_params_kwargs,
    "xlabel_kwargs": xlabel_kwargs,
    "ylabel_kwargs": ylabel_kwargs,
    "title_kwargs": title_kwargs,
    "legend_kwargs": legend_kwargs,
    "sim_kwargs": sim_kwargs,
    "subplot_kwargs": subplot_kwargs,
    "xlim": (transect_start, transect_end),
    "title": "Sacramento River between Sacramento and Yuba City",
    "ylim": (-75, 75),
    "xticks_kwargs": xticks_kwargs,
    "borelog_width": 1000,
    "streambed_bottom_col":"bathymetry_ft",
    "streambed_kwargs": streambed_kwargs
}

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
