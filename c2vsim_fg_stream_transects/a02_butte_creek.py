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
model_dir = os.path.join(data_dir, "c2vsimfg_v1.5")
preprocessor_dir = os.path.join(model_dir, "Preprocessor")
results_dir = os.path.join(model_dir, "Results")
aem_dir = os.path.join(data_dir, "sa7_supportingdata_20230428")

gwalloutfl_path = os.path.join(results_dir, "C2VSimFG_GW_HeadAll.out")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Preparing Stream and Groundwater Data from C2VSimFG v1.5
# In this section, we prepare the model data that will be used to make the stream transects. The workflow used hereafter
# consists of the following steps:
# 1. We load the Stream Specification File and the Nodal X-Y Coordinate File and generate a lines geodataframe
# of the streams represented in the model.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]

# Let's create the geodataframe of streams

# Path to stream definition file
streams_def_file_path = os.path.join(preprocessor_dir, "C2VSimFG_StreamsSpec.dat")

# Path to node definition file
nodes_def_file_path = os.path.join(preprocessor_dir, "C2VSimFG_Nodes.dat")

stream_lines_shp_path = os.path.join(data_dir, "stream_lines.shp")

stream_nodes_df = get_stream_nodes(streams_def_file_path, rating_points=10)

# We'll have to import the node definition file
gw_nodes_df = get_groundwater_nodes(nodes_def_file_path)

make_stream_lines_shp = True

if make_stream_lines_shp:

    streams_gdf = make_stream_lines_gdf(stream_nodes_df, gw_nodes_df)
    streams_gdf.to_file(stream_lines_shp_path)

else:
    streams_gdf = gpd.read_file(stream_lines_shp_path)

butte_creek_line_gdf = streams_gdf.loc[streams_gdf["name"].str.find("BUTTE")>=0].reset_index(drop=True)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# 2. Using the stream lines shapefile generated in the previous step, we calculate the distance of each stream node with
# respect to the origin of the stream line.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
# Let's add stream coordinates to stream nodes
stream_nodes_df = pd.merge(stream_nodes_df, gw_nodes_df, on="igw", how="left")

records = []
stream_line = list(streams_gdf.itertuples())[0]
for stream_line in streams_gdf.itertuples():
    stream_line_geometry = getattr(stream_line, "geometry")
    stream_name = getattr(stream_line, "name")
    for stream_node in stream_nodes_df[stream_nodes_df["name"] == stream_name].reset_index(drop=True).itertuples():
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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# 3. We load the Stratigraphy File and join ground surface elevations and layer thicknesses to the stream nodes.
# 4. We load the timeseries "Groundwater Head at All Nodes" file (C2VSimFG_GW_HeadAll.out) and join simulated groundwater heads
# to stream nodes.
# 5. For each node, we extract the groundwater heads for the highest active layer. To do so, we iterate through the
# layers and select the highest one that presents a thickness greater than zero and is above the bottom of the layer.
# 6. Select CASGEM wells and lithology logs within 500 m of Butte Creek.
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
# Let's load stratigraphy file into dataframe
stratigraphy_file_path = os.path.join(preprocessor_dir, "C2VSimFG_Stratigraphy.dat")


stratigraphy_df = get_stratigraphy(stratigraphy_file_path)

wt_stream_nodes_csv_path = os.path.join(results_dir, "gwallout_stream_nodes.csv")

make_wt_stream_nodes = True
if make_wt_stream_nodes:

    gwalloutfl_path: str = gwalloutfl_path
    date_width: int = 21
    head_width: int = 12
    header_lines: int = 5
    gwallout_df = load_gwalloutfl(gwalloutfl_path)
    # Let's convert heads from m to ft

    gwallout_df_long = pd.melt(gwallout_df, id_vars=["date", "layer"], var_name="igw", value_name="head_ft")

    # Let's add bottom elevations of layers
    stratigraphy_df["botm_lay_1"] = stratigraphy_df["gse"] - stratigraphy_df["thck_lay_1"]

    nlay = 4

    for lay in range(2, nlay+1):
        stratigraphy_df[f"botm_lay_{lay}"] = stratigraphy_df[f"botm_lay_{lay-1}"] - stratigraphy_df[f"thck_lay_{lay}"]

    gwallout_df_wide = pd.pivot(gwallout_df_long,index=["date", "igw"], columns=["layer"], values=["head_ft"])

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

    gwallout_stream_nodes_df["wte_ft"] = gwallout_stream_nodes_df.apply(get_wte, axis=1)

    gwallout_stream_nodes_df = gwallout_stream_nodes_df[
        ["date", "igw", "wte_ft", "head_ft_1", "head_ft_2", "head_ft_2", "head_ft_3", "head_ft_4"]
    ]

    gwallout_stream_nodes_df.to_csv(wt_stream_nodes_csv_path)
else:
    gwallout_stream_nodes_df = pd.read_csv(wt_stream_nodes_csv_path)

# %% [markdown] editable=true slideshow={"slide_type": ""}
# 6. Select CASGEM wells and AEM lithology logs within 500 m of Butte Creek
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
butte_creek_buffer = butte_creek_line_gdf.buffer(500)[0]
obs_wells_butte_creek_shp_path = os.path.join(data_dir, "obs_wells_butte_creek.shp")
lith_logs_butte_creek_shp_path = os.path.join(data_dir, "lith_logs_butte_creek.shp")
obs_butte_creek_csv_path = os.path.join(data_dir, "obs_wells_butte_creek.csv")

select_casgem_wells_and_lithology_logs = False
if select_casgem_wells_and_lithology_logs:
    # We download CASGEM wells for the counties that Butte Creek crosses
    url = "https://data.cnra.ca.gov/api/3/action/datastore_search_sql?"
    records = []
    for county in ['Butte', 'Sutter', 'Glenn', 'Tehama']:
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

    obs_wells_butte_creek_gdf = obs_wells_gdf.loc[obs_wells_gdf.geometry.within(butte_creek_buffer)].reset_index(drop=True)

    # Let's project the wells to the Butte Creek line
    obs_wells_butte_creek_gdf = project_points_gdf_to_line_string(
        obs_wells_butte_creek_gdf,
        butte_creek_line_gdf.loc[0,"geometry"],
    "site_code"
    )

    obs_wells_butte_creek_gdf.to_file(obs_wells_butte_creek_shp_path)

    lith_logs_shp_path = os.path.join(aem_dir, "WO7_HQ_LithologyWells.shp")

    lith_logs_gdf = gpd.read_file(lith_logs_shp_path)
    lith_logs_gdf = lith_logs_gdf.to_crs(streams_gdf.crs)

    lith_logs_butte_creek_gdf = lith_logs_gdf.loc[lith_logs_gdf.geometry.within(butte_creek_buffer)].reset_index(drop=True)

    # Let's project the lithology logs to the Butte Creek line
    lith_logs_butte_creek_gdf = project_points_gdf_to_line_string(
        lith_logs_butte_creek_gdf,
        butte_creek_line_gdf.loc[0,"geometry"],
    "WELLINFOID"
    )
    lith_logs_butte_creek_gdf.to_file(lith_logs_butte_creek_shp_path)

    # %% [markdown] editable=true slideshow={"slide_type": ""}
    # 7. Download CASGEM water level observations for the selected wells.
    # %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]

    label_col = "site_code"

    url = "https://data.cnra.ca.gov/api/3/action/datastore_search_sql?"
    dataset_code = "bfa9f262-24a1-45bd-8dc8-138bc8107266"
    site_codes_list = obs_wells_butte_creek_gdf[label_col].to_list()
    sql_query = f'''sql=SELECT * from "{
    dataset_code
    }" WHERE "site_code" IN ('{ "','".join(site_codes_list) }')'''
    sql_query = sql_query.replace(" ", "%20")
    sql_query = sql_query.replace('"', "%22")
    response_dmc = requests.get(url, params=sql_query)
    response_dict = json.loads(response_dmc.text)

    obs_butte_creek_df = pd.json_normalize(response_dict["result"]["records"])
    obs_butte_creek_df.to_csv(obs_butte_creek_csv_path, index=False)
else:
    obs_wells_butte_creek_gdf = gpd.read_file(obs_wells_butte_creek_shp_path)
    lith_logs_butte_creek_gdf = gpd.read_file(lith_logs_butte_creek_shp_path)
    obs_butte_creek_df = pd.read_csv(obs_butte_creek_csv_path)