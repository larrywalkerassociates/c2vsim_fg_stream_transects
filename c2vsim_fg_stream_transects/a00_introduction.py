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

from myst_nb import glue
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from shapely import LineString, Point
import shapely
import geopandas as gpd
import pandas as pd

from iwfm_lwa import (skip_until_flag, get_groundwater_nodes, load_gwalloutfl, get_stratigraphy, get_stream_nodes,
                      get_var, line_to_list)

warnings.filterwarnings('ignore')


# Path to the root of the repository
repo_dir = os.getcwd()
current_dir = os.getcwd()

while os.path.split(os.path.split(repo_dir)[0])[1] == "c2vsim_fg_stream_transects":
    repo_dir = os.path.split(repo_dir)[0]

os.chdir(repo_dir)

from functions.make_stream_lines_gdf import make_stream_lines_gdf
from functions.project_points_gdf_to_line_string import project_points_gdf_to_line_string

os.chdir(current_dir)

# Let's set key paths
data_dir = os.path.join(repo_dir, "data")
model_dir = os.path.join(data_dir, "c2vsimfg_v1.5")
preprocessor_dir = os.path.join(model_dir, "Preprocessor")
results_dir = os.path.join(model_dir, "Results")
aem_dir = os.path.join(data_dir, "sa7_supportingdata_20230428")

gwalloutfl_path = os.path.join(results_dir, "C2VSimFG_GW_HeadAll.out")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Methods
# Our approach consisted onIn this section, we prepare the model data that will be used to make the stream transects. The workflow used hereafter
# consists of the following steps:
# 1. We load the Stream Specification File and the Nodal X-Y Coordinate File and generate a lines geodataframe
# of the streams represented in the model.
# 2. Using the stream lines shapefile generated in the previous step, we calculate the distance of each stream node with
# respect to the origin of the stream line.
# 3. We load the Stratigraphy File and join ground surface elevations and layer thicknesses to the stream nodes.
# 4. We load the timeseries "Groundwater Head at All Nodes" file (C2VSimFG_GW_HeadAll.out) and join simulated groundwater heads
# to stream nodes.
# 5. For each node, we extract the groundwater heads for the highest active layer. To do so, we iterate through the
# layers and select the highest one that presents a thickness greater than zero and is above the bottom of the layer.
# 6. Select CASGEM wells and lithology logs within 500 m of Butte Creek.
# 6. Select CASGEM wells and AEM lithology logs within 500 m of Butte Creek
# 7. Download CASGEM water level observations for the selected wells.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
# Let's add stream coordinates to stream nodes
# Let's create the geodataframe of streams

# Path to stream definition file
streams_def_file_path = os.path.join(preprocessor_dir, "C2VSimFG_StreamsSpec.dat")

# Path to node definition file
nodes_def_file_path = os.path.join(preprocessor_dir, "C2VSimFG_Nodes.dat")

stream_lines_shp_path = os.path.join(data_dir, "stream_lines.shp")

#Buffer Distance in meters
buffer_distance_m = 500

stream_nodes_df = get_stream_nodes(streams_def_file_path, rating_points=10)

# We'll have to import the node definition file
gw_nodes_df = get_groundwater_nodes(nodes_def_file_path)

make_stream_lines_shp = False

if make_stream_lines_shp:

    streams_gdf = make_stream_lines_gdf(stream_nodes_df, gw_nodes_df)
    streams_gdf.to_file(stream_lines_shp_path)

else:
    streams_gdf = gpd.read_file(stream_lines_shp_path)

# Stream we will look at
stream_name = "SACRAMENTO"

stream_line_gdf = streams_gdf.loc[streams_gdf["name"]==stream_name].reset_index(drop=True)

stream_line_gdf = stream_line_gdf.drop(columns="id")

stream_line_gdf = stream_line_gdf.dissolve(by="name").reset_index()
stream_nodes_df = pd.merge(stream_nodes_df, gw_nodes_df, on="igw", how="left")

records = []

stream_line_geometry = stream_line_gdf.loc[0,"geometry"]

stream_nodes_df = stream_nodes_df[stream_nodes_df["name"] == stream_name].reset_index(drop=True)

for stream_node in stream_nodes_df.itertuples():
    stream_node_x = getattr(stream_node, "x")
    stream_node_y = getattr(stream_node, "y")
    irv = getattr(stream_node, "irv")
    point = Point(stream_node_x, stream_node_y)
    proj = stream_line_geometry.project(point)
    record = {
        "irv": irv,
        "proj": proj
    }
    records.append(record)

projs_df = pd.DataFrame(records)
stream_nodes_df = pd.merge(stream_nodes_df, projs_df, how="left")

def add_point_geometry(x):
    geometry = Point(x["x"], x["y"])
    return geometry

stream_nodes_df["geometry"] = stream_nodes_df.apply(add_point_geometry, axis=1)

save_stream_nodes = False
stream_nodes_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_stream_nodes.csv")
stream_nodes_shp_path = os.path.join(data_dir, f"{stream_name.lower()}_stream_nodes.shp")
if save_stream_nodes:
    stream_nodes_gdf = gpd.GeoDataFrame(
        stream_nodes_df,
        geometry=stream_nodes_df.geometry,
        crs=26910
    )
    stream_nodes_df.to_csv(stream_nodes_csv_path, index=False)
    stream_nodes_gdf.to_file(stream_nodes_shp_path)
else:
    stream_nodes_df = pd.read_csv(stream_nodes_csv_path)
    stream_nodes_gdf = gpd.read_file(stream_nodes_shp_path)


# Let's load stratigraphy file into dataframe
stratigraphy_file_path = os.path.join(preprocessor_dir, "C2VSimFG_Stratigraphy.dat")


stratigraphy_df = get_stratigraphy(stratigraphy_file_path)

wt_stream_nodes_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_gwallout_stream_nodes.csv")

make_wt_stream_nodes = False
if make_wt_stream_nodes:

    gwallout_df = load_gwalloutfl(gwalloutfl_path)
    # Let's convert heads from m to ft

    gwallout_df_long = pd.melt(gwallout_df, id_vars=["date", "layer"], var_name="igw", value_name="head_ft")

    gwallout_df = None

    # Let's add bottom elevations of layers
    stratigraphy_df["botm_lay_1"] = stratigraphy_df["gse"] - stratigraphy_df["thck_lay_1"]

    nlay = 4

    for lay in range(2, nlay+1):
        stratigraphy_df[f"botm_lay_{lay}"] = stratigraphy_df[f"botm_lay_{lay-1}"] - stratigraphy_df[f"thck_lay_{lay}"]

    gwallout_df_wide = pd.pivot(gwallout_df_long,index=["date", "igw"], columns=["layer"], values=["head_ft"])

    gwallout_df_long = None

    gwallout_df_wide = gwallout_df_wide.reset_index()

    colnames = gwallout_df_wide.columns

    colnames = [name[0]+"_"+str(name[1]) if len(str(name[1]))>0 else name[0] for name in colnames ]

    gwallout_df_wide.columns = colnames

    gwallout_df_wide = pd.merge(
        gwallout_df_wide,
    stratigraphy_df[
        [
            "igw", "gse"
        ] + [
            col for col in stratigraphy_df.columns if "botm" in col
        ] + [
            col for col in stratigraphy_df.columns if "thck_lay" in col
        ]
    ],
        how="left"
    )

    def get_wte(x):
        gse = x["gse"]
        wte = x["head_ft_1"]
        lay = 1
        botm = x["botm_lay_1"]
        thck = x["thck_lay_1"]
        while ((
                wte >= gse
        ) | (
                wte <= botm
        ) | (
                thck == 0
        ))& (
            lay < nlay
        ):
            lay = lay+1
            wte = x[f"head_ft_{lay}"]
            botm = x[f"botm_lay_{lay}"]
            thck = x[f"thck_lay_{lay}"]
        return wte

    gwallout_stream_nodes_df = pd.merge(
        gwallout_df_wide,
        stream_nodes_df["igw"],
        how="right"
    ).reset_index(drop=True)

    gwallout_df_wide = None

    gwallout_stream_nodes_df["wte_ft"] = gwallout_stream_nodes_df.apply(get_wte, axis=1)

    gwallout_stream_nodes_df = gwallout_stream_nodes_df[
        ["date", "igw", "wte_ft", "head_ft_1", "head_ft_2", "head_ft_2", "head_ft_3", "head_ft_4"]
    ]
    gwallout_stream_nodes_df = pd.merge(
        gwallout_stream_nodes_df,
        stream_nodes_df["igw"],

        how="right"
    ).reset_index(drop=True)

    gwallout_stream_nodes_df.to_csv(wt_stream_nodes_csv_path, index=False)
else:
    gwallout_stream_nodes_df = pd.read_csv(wt_stream_nodes_csv_path)


all_ts_df = pd.DataFrame(
    {"date": pd.to_datetime(gwallout_stream_nodes_df["date"].unique())}

)

all_ts_df["date_sim"] = all_ts_df["date"]
buffer = stream_line_gdf.buffer(buffer_distance_m)[0]
obs_wells_shp_path = os.path.join(data_dir, f"obs_wells_{stream_name.lower()}.shp")
lith_logs_shp_path = os.path.join(data_dir, f"lith_logs_{stream_name.lower()}.shp")
obs_csv_path = os.path.join(data_dir, f"obs_wells_{stream_name.lower()}.csv")

select_casgem_wells_and_lithology_logs = False
if select_casgem_wells_and_lithology_logs:
    # We download CASGEM wells for the counties that Butte Creek crosses
    url = "https://data.cnra.ca.gov/api/3/action/datastore_search_sql?"
    records = []
    for county in ['Butte', 'Colusa', 'Glenn', 'Sacramento', 'Shasta','Sutter', 'Tehama', 'Yolo']:
        sql_county = f'''sql=SELECT * from "af157380-fb42-4abf-b72a-6f9f98868077" WHERE "county_name" IN ('{county}')'''
        sql_county = sql_county.replace(" ", "%20")
        sql_county = sql_county.replace('"', "%22")
        response_dmc = requests.get(url, params=sql_county)

        response_dict = json.loads(response_dmc.text)
        records += response_dict["result"]["records"]

    obs_wells_df = pd.json_normalize(records)
    def add_geom_points(x):
        geometry = Point(x["longitude"], x["latitude"])
        return geometry

    obs_wells_df["geometry"] = obs_wells_df.apply(add_geom_points, axis=1)
    obs_wells_gdf = gpd.GeoDataFrame(obs_wells_df, geometry=obs_wells_df.geometry, crs=4326)

    obs_wells_gdf = obs_wells_gdf.to_crs(streams_gdf.crs)

    obs_wells_gdf = obs_wells_gdf.loc[obs_wells_gdf.geometry.within(buffer)].reset_index(drop=True)

    # Let's project the wells to the Sacramento River line
    obs_wells_gdf = project_points_gdf_to_line_string(
        obs_wells_gdf,
        stream_line_gdf.loc[0,"geometry"],
    "site_code"
    )



    lith_logs_shp_path_in = os.path.join(aem_dir, "WO7_HQ_LithologyWells.shp")

    lith_logs_gdf = gpd.read_file(lith_logs_shp_path_in)
    lith_logs_gdf = lith_logs_gdf.to_crs(streams_gdf.crs)

    lith_logs_gdf = lith_logs_gdf.loc[lith_logs_gdf.geometry.within(buffer)].reset_index(drop=True)

    # Let's project the lithology logs to the Butte Creek line
    lith_logs_gdf = project_points_gdf_to_line_string(
        lith_logs_gdf,
        stream_line_gdf.loc[0,"geometry"],
    "WELLINFOID"
    )
    lith_logs_gdf.to_file(lith_logs_shp_path)

    label_col = "site_code"

    url = "https://data.cnra.ca.gov/api/3/action/datastore_search_sql?"
    dataset_code = "bfa9f262-24a1-45bd-8dc8-138bc8107266"
    site_codes_list = obs_wells_gdf[label_col].to_list()
    records = []
    for site_code in site_codes_list:
        sql_query = f'''sql=SELECT * from "{
        dataset_code
        }" WHERE "site_code" IN ('{site_code}')'''
        sql_query = sql_query.replace(" ", "%20")
        sql_query = sql_query.replace('"', "%22")
        response_dmc = requests.get(url, params=sql_query)
        response_dict = json.loads(response_dmc.text)
        records += response_dict["result"]["records"]

    obs_df = pd.json_normalize(records)
    obs_df["date"] = pd.to_datetime(obs_df["msmt_date"])
    obs_df = obs_df.sort_values(by="date").reset_index(drop=True)
    obs_df = pd.merge_asof(
            obs_df, all_ts_df, on="date", direction="nearest"
        )
    obs_df = obs_df.drop(columns=["date"])
    obs_df = obs_df.rename(columns={
        "date_sim": "date"})
    # Let's average by well and date
    obs_df["gwe"] = obs_df["gwe"].astype(float)
    obs_df = obs_df[[label_col,"date","gwe"]].groupby(
        [label_col, "date"]
    ).mean().reset_index()
    obs_df = obs_df.loc[
        (
                obs_df["date"] >= pd.to_datetime("1973-09-01")
        ) & (obs_df["date"] <= pd.to_datetime("2021-09-30"))
        ].reset_index(drop=True)
    obs_df.to_csv(obs_csv_path, index=False)

    obs_wells_gdf = obs_wells_gdf.loc[
        obs_wells_gdf["site_code"].isin(obs_df["site_code"])
    ].reset_index(drop=True)

    obs_wells_gdf.to_file(obs_wells_shp_path)
else:
    obs_wells_gdf = gpd.read_file(obs_wells_shp_path)
    lith_logs_gdf = gpd.read_file(lith_logs_shp_path)
    obs_df = pd.read_csv(obs_csv_path)



# Let's get lithologies now
lith_csv_in_path = os.path.join(aem_dir, "AEM_WELL_LITHOLOGY_csv_WO7_20230327_HQonly.csv")
lith_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_lithology.csv")

make_lithology = False
if make_lithology:
    lith_df = pd.read_csv(lith_csv_in_path)

    lith_df = lith_df.rename(columns={"WELL_INFO_ID": "WELLINFOID"})

    lith_df["GROUND_SURFACE_ELEVATION_ft"] = lith_df["GROUND_SURFACE_ELEVATION_m"] * 3.28084
    lith_df["LITH_TOP_DEPTH_ft"] = lith_df["LITH_TOP_DEPTH_m"] * 3.28084
    lith_df["LITH_BOT_DEPTH_ft"] = lith_df["LITH_BOT_DEPTH_m"] * 3.28084

    # Let's select only the lithologies for Butte Creek
    lith_df = pd.merge(
        lith_df,
        lith_logs_gdf["WELLINFOID"],
        how="right",
    ).reset_index(drop=True)
    lith_df.to_csv(lith_csv_path, index=False)
else:
    lith_df = pd.read_csv(lith_csv_path)

# We get the axis limits so that all the points in the cross section are shown in the map

stream_line_string = stream_line_gdf.loc[0,"geometry"]

# Let's get the limits of the map
xsec_start = stream_line_string.length
xsec_end = 0

# Coordinates associated with the start and end of the cross section
xsec_start_coords = list(stream_line_string.interpolate(xsec_start).coords)[0]
xsec_end_coords = list(stream_line_string.interpolate(xsec_end).coords)[0]

x_min = xsec_start_coords[0]
x_max = xsec_end_coords[0]
y_min = xsec_start_coords[1]
y_max = xsec_end_coords[1]


for row in stream_nodes_gdf.itertuples():
    geometry = getattr(row, "geometry")
    point = geometry.coords[0]
    proj_d = stream_line_string.project(geometry)
    # We only consider points within the start and end of the cross section
    if (proj_d <= xsec_start) & (proj_d >= xsec_end):
        if point[0] < x_min:
            x_min = point[0]
        if point[0] > x_max:
            x_max = point[0]
        if point[1] < y_min:
            y_min = point[1]
        if point[1] > y_max:
            y_max = point[1]

# We set an offset of 500 m for axes limits
axis_lims_offset = 5000

step_miles_x = np.power(
    10, np.modf
            (np.log10(
            (
                    x_max - x_min
            ) / 1609.34
        )
        )[1]
)

xticks = np.arange(
            x_min-axis_lims_offset,
            x_max+axis_lims_offset,
            step_miles_x * 1609.34)

xticks_gdf = gpd.GeoDataFrame(
    {"utm_10N": xticks},
    geometry=[Point(x, y_min) for x in xticks],
    crs=26910
)

xticks_gdf = xticks_gdf.to_crs(4326)

# Now, we format the geographic coordinates for the tick labels into mins degrees, etc

def utm10n_to_geographic_lon(x):
    point = x["geometry"]
    lon = point.x
    if lon < 0:
        lon_deg = int(np.ceil(lon))
        lon_min = int(np.ceil((lon-lon_deg)*60))
    else:
        lon_deg = int(np.floor(lon))
        lon_min = int(np.floor((lon - lon_deg) * 60))
    lon_sec = (lon - lon_deg - lon_min / 60) * 3600
    if lon_deg < 0:
        lon_hem = "W"
    else:
        lon_hem = "E"

    lon_str = f"{abs(lon_deg)}°{abs(lon_min)}'{abs(lon_sec):.1f}\"{lon_hem}"
    return lon_str

xticks_gdf["lon"] = xticks_gdf.apply(utm10n_to_geographic_lon, axis=1)

step_miles_y = np.power(
    5, np.modf
            (np.log10(
            (
                    y_max - y_min
            ) / 1609.34
        )
        )[1]
)

yticks = np.arange(
            y_min,
            y_max,
            step_miles_y * 1609.34)

yticks_gdf = gpd.GeoDataFrame(
    {"utm_10N": yticks},
    geometry=[Point(x_min, y) for y in yticks],
    crs=26910
)

yticks_gdf = yticks_gdf.to_crs(4326)

def utm10n_to_geographic_lat(x):
    point = x["geometry"]
    lat = point.y
    if lat < 0:
        lat_deg = int(np.ceil(lat))
        lat_min = int(np.ceil((lat-lat_deg)*60))
    else:
        lat_deg = int(np.floor(lat))
        lat_min = int(np.floor((lat - lat_deg) * 60))
    lat_sec = (lat - lat_deg - lat_min / 60) * 3600
    if lat_deg < 0:
        lat_hem = "S"
    else:
        lat_hem = "N"

    lat_str = f"{abs(lat_deg)}°{abs(lat_min)}'{abs(lat_sec):.1f}\"{lat_hem}"
    return lat_str


yticks_gdf["lat"] = yticks_gdf.apply(utm10n_to_geographic_lat, axis=1)

# We will create a geodataframe with the tickmarks, which we will project to EPSG:4326
# the tickmarks in geographic coordinates


subplot_kwargs = {"ncols": 2,
                  "width_ratios": [0.9, 0.1],
                  "figsize": (12, 8)}

fig, axes = plt.subplots(**subplot_kwargs)

# Let's plot the Sacramento River stream line now

stream_line_style_kwargs = {"lw": 2}
legend_elements = []

axes[0] = stream_line_gdf.plot(ax=axes[0], color="blue", **stream_line_style_kwargs)
legend_elements.append(
    Line2D([0], [0], color="blue", label= "Sacramento River", **stream_line_style_kwargs)
)

# Let's plot the CASGEM wells now
casgem_well_style_kwargs = {}
legend_elements.append(Line2D([0], [0], label='CASGEM Wells',
                              markerfacecolor="#CB6015",marker="o",
                              color='w', **casgem_well_style_kwargs))
axes[0] = obs_wells_gdf.plot(ax=axes[0], color = "#CB6015", **casgem_well_style_kwargs)

axes[1].set_axis_off()

# Let's plot the lithology logs now
lithology_log_style_kwargs = {"marker": "^"}
legend_elements.append(Line2D([0], [0], label='Lithology Logs',
                              markerfacecolor="#A2495E",
                              color='w', **lithology_log_style_kwargs))
axes[0] = lith_logs_gdf.plot(ax=axes[0], color = "#A2495E", **lithology_log_style_kwargs)

axes[1].legend(handles=legend_elements)

axes[0].set_xlim(x_min-axis_lims_offset, x_max+axis_lims_offset)
axes[0].set_ylim(y_min-axis_lims_offset, y_max+axis_lims_offset)

axes[0].set_xticks(xticks, labels=xticks_gdf["lon"].to_list(), rotation=90)
axes[0].set_yticks(yticks, labels=yticks_gdf["lat"].to_list())

plt.tight_layout()

#glue("scenario_name", scenario_name, display=False)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} f_dmc039
# :figwidth: 800px
# :name: "fig:f_dmc039"
#
# Simulated NVRRWP fraction at node DMC039 for the {glue:}`scenario_name` scenario.
# ```
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
