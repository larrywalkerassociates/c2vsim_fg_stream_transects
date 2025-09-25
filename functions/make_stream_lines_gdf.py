""" Function to create a GeoDataframe of stream lines from a DataFrame of stream nodes and another DataFrame of
groundwater nodes created with functions get_stream_nodes and get_groundwater_nodes from iwfm_lwa package"""


from shapely import LineString
import geopandas as gpd
import pandas as pd

def make_stream_lines_gdf(
        stream_nodes_df: pd.DataFrame,
        gw_nodes_df: pd.DataFrame,
        igw_col: str="igw",
        id_col: str= "id",
        x_col: str="x",
        y_col: str="y",
        name_col: str="name",
        crs: int=26910,
)-> gpd.GeoDataFrame:

    """
    Create GeoDataframe of stream lines from a DataFrame of stream nodes and another DataFrame of groundwater nodes

    This function creates a GeoDataframe of stream lines from a DataFrame of stream nodes and another DataFrame of
    groundwater nodes created with functions get_stream_nodes and get_groundwater_nodes from iwfm_lwa package.

    Parameters
    ----------
    stream_nodes_df : pd.DataFrame
        Dataframe of stream nodes created using functions get_stream_nodes from iwfm_lwa package.
    gw_nodes_df : pd.DataFrame
        Dataframe of groundwater nodes created using functions get_groundwater_nodes from iwfm_lwa package.
    igw_col : str, optional
        Column containing groundwater node index for stream_nodes_df and gw_nodes_df.
        Default is "igw".
    id_col : str, optional
        Column containing stream ids for stream_nodes_df.
        Default is "id".
    x_col: str, optional
        Column containing node x coordinates for gw_nodes_df.
        Default is "x".
    y_col: str, optional
        Column containing node y coordinates for gw_nodes_df.
        Default is "y".
    crs: int, optional
        EPSG code of reference coordinate system used in the model.
        Default is 6414 (NAD83 / California Albers in m).

    Returns
    -------
    geopandas.GeoDataFrame
        A table with one row per stream and the following columns:
        - 'id' (int): Stream ID.
        - 'name' (str): Stream name.

    """

    # We join coordinates to stream nodes
    stream_nodes_df = pd.merge(stream_nodes_df, gw_nodes_df, on=igw_col, how="left")

    def add_coords(x):
        coords = (x[x_col], x[y_col])
        return coords

    stream_nodes_df["coords"] = stream_nodes_df.apply(add_coords, axis=1)

    # Let's make line with list of coordinates for each node
    reaches = stream_nodes_df[id_col].unique()

    records = []

    for reach in reaches:
        stream_network_df_local = stream_nodes_df.loc[
            stream_nodes_df[id_col] == reach
            ].reset_index(drop=True)
        geometry = LineString(stream_network_df_local.coords.to_list())

        names = stream_network_df_local[name_col].unique()
        if len(names) > 1:
            print(f"Warning: Multiple names found for stream ID {reach}: {names}")

        record = {
            id_col: reach,
            name_col: names[0],
            "geometry": geometry,
        }

        records.append(record)

    streams_df = pd.DataFrame(records)

    streams_gdf = gpd.GeoDataFrame(
        data=streams_df, geometry="geometry", crs=crs)

    return streams_gdf