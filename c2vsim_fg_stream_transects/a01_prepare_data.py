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
import rasterio
import requests
import warnings

from myst_nb import glue
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib_map_utils.core.north_arrow import NorthArrow, north_arrow
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
bathymetry_path = os.path.join(data_dir, "sacramento_river_bathymetry_2m.tif")
sacramento_path = os.path.join(data_dir, "sacramento.shp")
yuba_city_path = os.path.join(data_dir, "yuba_city.shp")
chico_path = os.path.join(data_dir, "chico.shp")
red_bluff_path = os.path.join(data_dir, "red_bluff.shp")

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Methods
# We extracted stream nodes and their corresponding simulated groundwater heads from CV2SIM-FG v1.5
# {cite}`california_department_of_water_resources_california_2025-1`,observed groundwater levels from CASGEM
# {cite}`california_department_of_water_resources_california_2025`, and high-quality lithology logs from DWR's Airborne Electromagnetic Survey (AEM) for areas 6
# {cite}`california_department_of_water_resources_california_2023` and 7
# {cite}`california_department_of_water_resources_california_2023-1`
#
# We connected the stream nodes to create a Sacramento River streamline. For each node and timestep, we computed
# groundwater table as the head of the shallowest layer where the head exceeded the bottom elevation. Using simulated
# heads and CV2SIM-FG v1.5 bottom layer elevations, we calculated the groundwater table elevation for each node and timestep as the head of the highest layer where the head is above
# the layer bottom elevation. We aggregated AEM lithology textures into broader categories for visualization, as shown in
# {numref}`lith_plot_cats`. We included CASGEM and AEM lithology wells within 500 m (1,640 ft) of the river. Where CASGEM
# well labels overlapped, we prioritized observation over irrigation and domestic wells, and shallow over deep wells.
# When shallow observation wells overlapped, we manually adjusted the x-coordinate of the projection onto the transect
# for label readability.
#
# ```{table} Plotting categories for AEM lithology logs.
# :name: lith_plot_cats
# | Texture                 | Plotting Group |
# | ----------------------- | -------------- |
# | Boulders                | Gravel         |
# | Clay                    | Clay           |
# | Cobbles                 | Gravel         |
# | Gravel                  | Gravel         |
# | Rock - Intrusive        | Rock           |
# | Rock - Metamorphic      | Rock           |
# | Rock - Sedimentary      | Rock           |
# | Rock - Undifferentiated | Rock           |
# | Rock - Volcanic         | Rock           |
# | Sand                    | Sand           |
# | Silt                    | Silt           |
# | Soil                    | Soil           |
# | Unknown                 | Unknown        |
# ```
#
# We analyzed a 190-mile section of the Sacramento River from Red Bluff to Sacramento, which we divided in three sections
# for clarity. We extracted Sacramento River bathymetry at stream nodes from {cite}`yurok_tribe_topo-bathymetric_2024`
# and {cite}`california_department_of_water_resources_revised_2019`. We assigned cross section mileage relative to the
# confluence with the San Joaquin River. The soutern transect covers from Sacramento to Yuba City (mile 50 to 120).
# The central transect spans from Yuba City to Chico (mile 120 to 190). The northern transect ranges from Chico to Red
# Bluff (mile 190 to 240).
#
# For our workflow, we used libraries rasterio {cite}`gillies_rasterio_2013`, shapely
# {cite}`gillies_shapely_2025`, pandas {cite}`mckinney_data_2010`, geopandas {cite}`jordahl_geopandas_2020`, and
# numpy {cite}`harris_array_2020` for data processing. For preparing the animations we used matplotlib
# {cite}`hunter_matplotlib_2007`. We handled deployment was handled using jupyter book
# {cite}`executable_books_community_jupyter_2020` and GitHub Pages {cite}`github_github_nodate`.
#
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
    values = []
    with rasterio.open(bathymetry_path) as ds:
        vals_array = ds.read(1)
        for point in stream_nodes_gdf["geometry"]:
            x = point.x
            y = point.y
            row,col = ds.index(x,y)
            if (row>=0) & (row<vals_array.shape[0]) & (col>=0) & (col<vals_array.shape[1]):
                value = vals_array[row,col]
            else:
                value = np.nan
            values.append(value)
    # let's remove 0s, which QGIS used to represent nas
    values = [val if val != 0 else np.nan for val in values]
    # Let's convert from m to ft
    values = [val * 3.28084 if not np.isnan(val) else np.nan for val in values]
    # Let's add bathymetry values to stream nodes
    stream_nodes_gdf["bathymetry_ft"] = values
    stream_nodes_df["bathymetry_ft"] = values

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

    lith_logs_shp_path_in_sa6 = os.path.join(aem_dir, "WO6_HQ_LithologyWells.shp")

    lith_logs_gdf = gpd.read_file(lith_logs_shp_path_in)
    lith_logs_gdf = lith_logs_gdf.to_crs(streams_gdf.crs)
    lith_logs_gdf["WELLINFOID"] = lith_logs_gdf["WELLINFOID"].astype(str)+"_sa7"
    lith_logs_sa6_gdf = gpd.read_file(lith_logs_shp_path_in_sa6)
    lith_logs_sa6_gdf = lith_logs_sa6_gdf.to_crs(streams_gdf.crs)
    lith_logs_sa6_gdf["WELLINFOID"] = lith_logs_sa6_gdf["WELLINFOID"].astype(str) + "_sa6"


    lith_logs_gdf = pd.concat(
        [
            lith_logs_gdf,
            lith_logs_sa6_gdf
        ],
        ignore_index=True
    )



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
lith_csv_in_path_6 = os.path.join(aem_dir, "AEM_WELL_LITHOLOGY_csv_WO6_20230103_HQonly.csv")
lith_csv_path = os.path.join(data_dir, f"{stream_name.lower()}_lithology.csv")

make_lithology = False
if make_lithology:
    lith_df = pd.read_csv(lith_csv_in_path)
    lith_df = lith_df.rename(columns={"WELL_INFO_ID": "WELLINFOID"})
    lith_df["WELLINFOID"] = lith_df["WELLINFOID"].astype(str) + "_sa7"

    lith_df_6 = pd.read_csv(lith_csv_in_path_6)
    lith_df_6 = lith_df_6.rename(columns={"WELL_INFO_ID": "WELLINFOID"})
    lith_df_6["WELLINFOID"] = lith_df_6["WELLINFOID"].astype(str) + "_sa6"


    lith_df = pd.concat(
        [
            lith_df,
            lith_df_6
        ],
        ignore_index=True
    )



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

xlims_utm10n_gdf = gpd.GeoDataFrame(
    {"utm_10N": [x_min-axis_lims_offset, x_max+axis_lims_offset]},
    geometry=[Point(x_min-axis_lims_offset, y_min), Point(x_max+axis_lims_offset, y_min)],
    crs=26910
)

ylims_utm10n_gdf = gpd.GeoDataFrame(
    {"utm_10N": [y_min-axis_lims_offset, y_max+axis_lims_offset]},
    geometry=[Point(x_min, y_min-axis_lims_offset), Point(x_min, y_max+axis_lims_offset)],
    crs=26910
)


# Let's convert to geographic coordinates
xlims_geographic_gdf = xlims_utm10n_gdf.to_crs(4326)
x_max_degrees = xlims_geographic_gdf.loc[1, "geometry"].x
x_min_degrees = xlims_geographic_gdf.loc[0, "geometry"].x

ylims_geographic_gdf = ylims_utm10n_gdf.to_crs(4326)
y_max_degrees = ylims_geographic_gdf.loc[1, "geometry"].y
y_min_degrees = ylims_geographic_gdf.loc[0, "geometry"].y

step_degrees_x = np.modf(x_max_degrees - x_min_degrees)[1]/4

step_degrees_y = np.modf(y_max_degrees - y_min_degrees)[1]/4

xticks = np.arange(
            np.modf(x_min_degrees/step_degrees_x)[1]*step_degrees_x,
            np.modf(x_max_degrees/step_degrees_x)[1]*step_degrees_x,
            step_degrees_x)

yticks = np.arange(
            np.modf(y_min_degrees/step_degrees_y)[1]*step_degrees_y,
            np.modf(y_max_degrees/step_degrees_y)[1]*step_degrees_y,
            step_degrees_y)

xticks_labels = [f"{np.abs(np.modf(tick)[1]):.0f}°{np.modf(np.modf(tick)[0]*60)[1]:.0f}'{np.modf(np.modf(tick)[0]*60)[0]*3600:.1f}\" W" for tick in xticks]

yticks_labels = [f"{np.modf(tick)[1]:.0f}°{np.modf(np.modf(tick)[0]*60)[1]:.0f}'{np.modf(np.modf(tick)[0]*60)[0]*3600:.1f}\" N" for tick in yticks]

step_miles_y = np.power(
    5, np.modf
            (np.log10(
            (
                    y_max - y_min
            ) / 1609.34
        )
        )[1]
)

xticks_gdf = gpd.GeoDataFrame(
    {"4326": xticks},
    geometry=[Point(x, xlims_geographic_gdf.loc[0, "geometry"].y) for x in xticks],
    crs=4326
)

yticks_gdf = gpd.GeoDataFrame(
    {"4326": yticks},
    geometry=[Point(ylims_geographic_gdf.loc[0, "geometry"].x, y) for y in yticks],
    crs=4326
)

xticks_gdf = xticks_gdf.to_crs(26910)

yticks_gdf = yticks_gdf.to_crs(26910)

xticks_gdf["x"] = xticks_gdf["geometry"].x

yticks_gdf["y"] = yticks_gdf["geometry"].y

# We will create a geodataframe with the tickmarks, which we will project to EPSG:4326
# the tickmarks in geographic coordinates

sacramento_gdf = gpd.read_file(sacramento_path)
yuba_city_gdf = gpd.read_file(yuba_city_path)
chico_gdf = gpd.read_file(chico_path)
red_bluff_gdf = gpd.read_file(red_bluff_path)

# Let's make geodataframe with transects postmiles
postmiles = [50, 100, 161, 231]

# Convert to distance from confluence with San Joaquin
postmiles_text = [int(np.modf((stream_line_string.length/1609.34-postmile)/10)[1]*10) for postmile in postmiles]

postmiles_pts = [stream_line_string.interpolate(postmile*1609.34) for postmile in postmiles]

postmiles_gdf = gpd.GeoDataFrame(
    data={"postmile": postmiles_text,
          "geometry": postmiles_pts},
    crs=26910
)

subplot_kwargs = {"ncols": 2,
                  "width_ratios": [0.9, 0.1],
                  "figsize": (6, 8)}

fig, axes = plt.subplots(**subplot_kwargs)

# Let's plot the Sacramento River stream line now

stream_line_style_kwargs = {"lw": 2}
legend_elements = []

axes[0] = stream_line_gdf.plot(ax=axes[0], color="blue", **stream_line_style_kwargs)
legend_elements.append(
    Line2D([0], [0], color="blue", label= "Sacramento River", zorder = 1, **stream_line_style_kwargs)
)

# Let's plot the CASGEM wells now
casgem_well_style_kwargs = {}
legend_elements.append(Line2D([0], [0], label='CASGEM Wells',
                              markerfacecolor="#CB6015",marker="o",
                              color='w', markersize=10, **casgem_well_style_kwargs))
axes[0] = obs_wells_gdf.plot(ax=axes[0], color = "#CB6015", zorder = 2, **casgem_well_style_kwargs)

axes[1].set_axis_off()

# Let's plot the lithology logs now
lithology_log_style_kwargs = {"marker": "^"}
axes[0] = lith_logs_gdf.plot(ax=axes[0], color = "#9BAA4B", zorder = 2, **lithology_log_style_kwargs)

legend_elements.append(Line2D([0], [0], label='Lithology Logs',
                              markerfacecolor="#9BAA4B",
                              color='w', markersize=10, **lithology_log_style_kwargs))

axes[0] = sacramento_gdf.plot(ax=axes[0], color = 'none',  hatch = "////", edgecolor="#555555")
legend_elements.append(Patch(facecolor="none", hatch = "////", edgecolor="#555555", label="Sacramento"))

axes[0] = yuba_city_gdf.plot(ax=axes[0], color = 'none',  hatch = "\\\\\\\\", edgecolor="#F6B704")
legend_elements.append(Patch(facecolor="none", hatch = "\\\\\\\\", edgecolor="#F6B704", label="Yuba City"))

axes[0] = chico_gdf.plot(ax=axes[0], color = 'none',  hatch = "----", edgecolor="#4B7164")
legend_elements.append(Patch(facecolor="none",  hatch = "----", edgecolor="#4B7164", label="Chico"))

axes[0] = red_bluff_gdf.plot(ax=axes[0], color = 'none',  hatch = "||||", edgecolor="#A2495E")
legend_elements.append(Patch(facecolor="none",  hatch = "||||", edgecolor="#A2495E", label="Red Bluff"))

axes[1].legend(handles=legend_elements)

axes[0].set_xlim(x_min-axis_lims_offset, x_max+axis_lims_offset)
axes[0].set_ylim(y_min-axis_lims_offset, y_max+axis_lims_offset)

axes[0].set_xticks(xticks_gdf["x"], labels=xticks_labels, rotation=90)
axes[0].set_yticks(yticks_gdf["y"], labels=yticks_labels)

asb = AnchoredSizeBar(axes[0].transData,
                      5*1609.34,
                      "5 mi",
                      size_vertical=1609.34,
                      loc="upper right",
                      pad=0.1,
                      borderpad=0.5,
                      sep=5,
                      frameon=False)
axes[0].add_artist(asb)
north_arrow(
    axes[0], location="lower left", rotation={"crs": stream_line_gdf.crs, "reference": "center"}
)

postmile_240 = postmiles_gdf.loc[0, "postmile"]
geom_pm_240 = postmiles_gdf.loc[0, "geometry"]
xy_coords_pm_240 = list(geom_pm_240.coords)[0]
axes[0].annotate(
        f"{postmile_240}",
        xy=xy_coords_pm_240,
        xytext=(xy_coords_pm_240[0]+8000, xy_coords_pm_240[1]),
        bbox=dict(boxstyle="square", fc="w"),
    )

postmile_190 = postmiles_gdf.loc[1, "postmile"]
geom_pm_190 = postmiles_gdf.loc[1, "geometry"]
xy_coords_pm_190 = list(geom_pm_190.coords)[0]
axes[0].annotate(
        f"{postmile_190}",
        xy=xy_coords_pm_190,
        xytext=(xy_coords_pm_190[0]-25000, xy_coords_pm_190[1]),
        bbox=dict(boxstyle="square", fc="w"),
    )

postmile_120 = postmiles_gdf.loc[2, "postmile"]
geom_pm_120 = postmiles_gdf.loc[2, "geometry"]
xy_coords_pm_120 = list(geom_pm_120.coords)[0]
axes[0].annotate(
        f"{postmile_120}",
        xy=xy_coords_pm_120,
        xytext=(xy_coords_pm_120[0]-25000, xy_coords_pm_120[1]),
        bbox=dict(boxstyle="square", fc="w"),
    )

postmile_60 = postmiles_gdf.loc[3, "postmile"]
geom_pm_60 = postmiles_gdf.loc[3, "geometry"]
xy_coords_pm_60 = list(geom_pm_60.coords)[0]
axes[0].annotate(
        f"{postmile_60}",
        xy=xy_coords_pm_60,
        xytext=(xy_coords_pm_60[0]-25000, xy_coords_pm_60[1]),
        bbox=dict(boxstyle="square", fc="w"),
    )

# for row in postmiles_gdf.itertuples():
#     postmile = getattr(row, "postmile")
#     geom = getattr(row, "geometry")
#     axes[0].annotate(
#         f"{postmile}",
#         xy=list(geom.coords)[0],
#         xytext=list(geom.coords)[0],
#         bbox=dict(boxstyle="square", fc="w"),
#     )

plt.tight_layout()

glue("sac_river_fig", fig, display=False)

plt.close()
# %% [markdown] editable=true slideshow={"slide_type": ""}
# ```{glue:figure} sac_river_fig
# :figwidth: 800px
# :name: "sac_river_fig"
#
# Sacramento River and selected CASGEM wells and lithology logs. Distance along the Sacramento River from the confluence
# with the San Joaquin River at transect sections boundaries are shown in white boxes in units of miles.
# ```
# %% editable=true slideshow={"slide_type": ""} tags=["remove-input"]
