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
# # Sacramento River: R
# In this section, we will make a stream transect for Butte Creek using the C2VSim-FG model streambed and groundwater
# table elevations, DWR CASGEM groundwater levels observations, and DWR's AEM lithology logs.
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

# We remove deep piezometers of nested well 398028N1220294
obs_lut_df = obs_lut_df[~obs_lut_df["site_code"].isin(
    ["398028N1220294W001", "398028N1220294W002", "398028N1220294W003"]
)].reset_index(drop=True)

# We remove deep well 401774N1222253W001 in crowded area
obs_lut_df = obs_lut_df[obs_lut_df["site_code"]!="401774N1222253W001"].reset_index(drop=True)

# We remove well without depth information 401678N1222236W001 in crowded area
obs_lut_df = obs_lut_df[obs_lut_df["site_code"]!="401678N1222236W001"].reset_index(drop=True)

# We remove too shallow well 401728N1222223W001 in crowded area
obs_lut_df = obs_lut_df[obs_lut_df["site_code"]!="401728N1222223W001"].reset_index(drop=True)

# We remove irrigation well 401550N1222006W001 in crowded area
obs_lut_df = obs_lut_df[obs_lut_df["site_code"]!="401550N1222006W001"].reset_index(drop=True)

# Let's load the observed groundwater levels

obs_df = pd.read_csv(obs_wells_csv_path, parse_dates=["date"])

obs_df.loc[(obs_df["site_code"]=="396568N1217818W001")&(obs_df["date"]=="2015-09-30"),"date"] = pd.to_datetime(
    "2015-10-31"
)


step_miles = 1 * np.power(
    10, np.modf
            (np.log10(
            (
                    stream_nodes_df["proj"].max() - stream_nodes_df["proj"].min()
            ) / 1609.34
        )
        )[1] - 1
)

transect_start =100*1609.34  # 100 miles from north end, Chico

transect_end = 50*1609.34 # 50 miles from north end, Red Bluff

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
            "color": "w"
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
    "title": "Sacramento River between Red Bluff and Chico",
    "ylim": (100, 250),
    "xticks_kwargs": xticks_kwargs,
    "borelog_width": 1000,
    "streambed_bottom_col":"bathymetry_ft",
    "streambed_kwargs": streambed_kwargs
}

# html_video = make_transect_animation(
#         stream_nodes_df,
#         sim_df,
#         obs_df,
#         obs_lut_df,
#         lith_lut_df,
#         lith_df,
#     **kwargs)

###################################################################################################################
date_col = kwargs.pop("date_col", "date")
obs_col = kwargs.get("obs_col", "gwe")
obs_kwargs = kwargs.get("obs_kwargs", None)
proj_d_col = kwargs.get("proj_d_col", "proj")
sim_col = kwargs.get("sim_col", "wte_ft")
sim_kwargs = kwargs.get("sim_kwargs", None)
title = kwargs.get("title", "Butte Creek Stream Transect")
title_kwargs = kwargs.get("title_kwargs", None)
well_id_col = kwargs.get("well_id_col", "site_code")
node_key_col = kwargs.pop("node_key_col", "igw")
interval = kwargs.pop("interval", 500)  # milliseconds between frames

if obs_kwargs is None:
    obs_kwargs = {
        "label": 'Observed GWE',
        "markerfacecolor": "#4B7164",
        "marker": 'o',
        "color": "w"
    }

if sim_kwargs is None:
    sim_kwargs = {
        "color": "navy",
        "label": "Simulated GWE"
    }
    kwargs["sim_kwargs"] = sim_kwargs

if title_kwargs is None:
    title_kwargs = {"pad": 150}
    kwargs["title_kwargs"] = title_kwargs

fig, axes = make_transect_layout(
    stream_nodes_df,
    sim_df,
    obs_df,
    obs_lut_df,
    lith_lut_df,
    lith_df, **kwargs)

sim_kwargs.pop("label")
obs_kwargs.pop("label")
obs_kwargs["color"] = obs_kwargs.pop("markerfacecolor")


# Now, we make the update function that will make the frames in the animation
def update(
        ts,
        obs_df,
        obs_lut_df,
        sim_df,
        date_col=date_col,
        proj_d_col=proj_d_col,
        sim_col=sim_col,
        title=title,
        title_kwargs=title_kwargs,
        obs_kwargs=obs_kwargs,
        well_id_col=well_id_col):
    artists = []
    # List of artists labels that we need to remove every timestep
    artists_to_remove = obs_lut_df[well_id_col].tolist() + [sim_col] + [title]
    # Let's remove artists from previous frame
    for artist in axes[0].get_children():
        artist_label = artist.get_label()
        if artist_label in artists_to_remove:
            try:
                artist.remove()
            except ValueError:
                pass

    # Let's add observations for the given timestep
    obs_df_ts = obs_df[obs_df[date_col] == ts].reset_index(drop=True)
    for well in obs_df_ts[well_id_col].unique().tolist():
        x = obs_lut_df.loc[obs_lut_df[well_id_col] == well, proj_d_col].values[0]
        obs = obs_df_ts.loc[obs_df_ts[well_id_col] == well, obs_col].values[0]
        artists.append(axes[0].scatter(x, obs, label=well, **obs_kwargs))

    # Let's add simulated groundwater levels for the given timestep
    sim_df_ts = sim_df[sim_df[date_col] == ts].reset_index(drop=True)
    # Let's join with stream nodes to get the projection distance
    sim_df_ts = pd.merge(
        sim_df_ts,
        stream_nodes_df[[proj_d_col, node_key_col]],
        on=node_key_col,
        how="left")
    artists.append(
        axes[0].plot(
            sim_df_ts[proj_d_col],
            sim_df_ts[sim_col],
            label=sim_col,
            **sim_kwargs
        )
    )

    # let's set the new title label for the given timestep
    title = f"{ts.strftime('%Y %B')}"
    artists_to_remove[len(artists_to_remove) - 1] = title
    artists.append(axes[0].set_title(title, **title_kwargs))

    return artists


ani = animation.FuncAnimation(
    fig=fig,
    func=functools.partial(
        update,
        obs_df=obs_df,
        obs_lut_df=obs_lut_df,
        sim_df=sim_df
    ),
    frames=sim_df[date_col].unique(),
    interval=interval,
)

html_video = ani.to_html5_video()

plt.close()
###################################################################################################################

widgets.HTML(
        value=html_video,
        placeholder="",
        description="",
    )
