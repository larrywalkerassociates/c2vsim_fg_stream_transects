""" Function to project the geometries of a GeoDataFrame of points to a shapely LineString object"""
import logging

import geopandas as gpd
import pandas as pd
import shapely

def project_points_gdf_to_line_string(
        points_gdf: gpd.GeoDataFrame,
        line_geometry: shapely.LineString,
        id_col: str,
)-> gpd.GeoDataFrame:

    """

Project point geometries of a GeoDataFrame to a LineString and compute their
    position along the line.

    Parameters
    ----------
    points_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing Point geometries.
    line_geometry : shapely.LineString
        The line onto which points will be projected.

    Returns
    -------
    geopandas.GeoDataFrame
        A copy of `points_gdf` with an added column:
        - 'proj' (float): Distance from the line's start (units of CRS)


    Notes
    -----
    - Distances are computed in the units of the CRS. If your CRS is geographic (e.g., EPSG:4326),
      distances will be in degrees. Reproject to a suitable projected CRS (meters/feet) before using.
    - Points outside the line's extents project to the nearest point on the line (ends included).

    Raises
    ------
    ValueError
        If points_gdf does not contain a geometry column, if geometries are not Points.
    """


    if "geometry" not in points_gdf:
        raise ValueError("points_gdf must have a 'geometry' column.")


    if not points_gdf.geometry.geom_type.isin(["Point"]).all():
        raise ValueError("points_gdf must contain only Point geometries.")


    if points_gdf.geometry.is_empty.any():
        logging.warning("Some point geometries are empty; their projections will be NaN.")


    # Warn if CRS might lead to confusing units
    if hasattr(points_gdf, "crs") and points_gdf.crs is not None:
        try:
            if getattr(points_gdf.crs, "is_geographic", False):
                logging.warning(
                    "points_gdf CRS is geographic (degrees). Consider converting to a projected CRS "
                    "(e.g., UTM) to get distances in meters.",
                )
        except Exception:
            # Some CRS objects may not have is_geographic; ignore if not available
            pass
    else:
        logging.warning("points_gdf has no CRS; distance units may be ambiguous.")



    records = []
    for point in points_gdf.itertuples():
        point_geom = getattr(point, "geometry")
        pid = getattr(point, id_col)
        proj = line_geometry.project(point_geom)
        record = {
            id_col: pid,
            "proj": proj
        }
        records.append(record)

    projs_df = pd.DataFrame(records)
    points_gdf = pd.merge(points_gdf, projs_df, how="left")

    return points_gdf