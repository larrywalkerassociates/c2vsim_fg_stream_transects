""" Function to make the static layout of a transect figure"""

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

from functions.add_bore_logs import add_bore_logs


def make_transect_layout(
        stream_nodes_df,
        sim_df,
        obs_df,
        obs_lut_df,
        lith_lut_df,
        lith_df,
        **kwargs
)-> tuple[matplotlib.axes.Axes | np.ndarray, matplotlib.figure.Figure]:

    """
    Makes the static layout of a transect figure containing borehole lithology logs and streambed bottom elevations.

    Parameters
    ----------
    stream_nodes_df: pandas.DataFrame
        Dataframe containing stream node information, including projected distances and streambed bottom elevations.
    sim_df: pandas.DataFrame
        Dataframe containing simulated groundwater elevations at the stream nodes.
    obs_df: pandas.DataFrame
        Dataframe containing observed groundwater elevations at monitoring wells.
    obs_lut_df: pandas.DataFrame
        Dataframe containing lookup information for monitoring wells, including projected distances along the transect.
    lith_lut_df: pandas.DataFrame
        Dataframe containing lookup information for borehole lithology logs, including projected distances along the
        transect.
    lith_df: pandas.DataFrame
        Dataframe containing lithology log data for boreholes.
    borelog_id_col: str, optional
        Column name for borelog IDs in lith_logs_gdf and lith_df.
        Default is "WELLINFOID".
    borelog_width: float, optional
        Width of the borelog patches in the plot.
        Default is 1000.
    legend_items: list, optional
        List of matplotlib Line2D objects to include in the legend.
        If None, a default legend is created.
        Default is None.
    legend_kwargs: dict, optional
        Additional keyword arguments for the legend.
        Default is None.
    obs_col: str, optional
        Column name for observed groundwater elevations in obs_df.
        Default is "gwe".
    proj_d_col: str, optional
        Column name for projected distances along the transect in stream_nodes_df, obs_lut_df, and lith_lut_df.
        Default is "proj".
    sim_col: str, optional
        Column name for simulated groundwater elevations in sim_df.
        Default is "wte_ft".
    sim_kwargs: dict, optional
        Additional keyword arguments for plotting simulated groundwater elevations.
        Default is None.
    streambed_bottom_col: str, optional
        Column name for streambed bottom elevations in stream_nodes_df.
        Default is "botr".
    streambed_kwargs: dict, optional
        Additional keyword arguments for plotting streambed bottom elevations.
        Default is None.
    subplot_kwargs: dict, optional
        Additional keyword arguments for creating subplots.
        Default is None.
    title: str, optional
        Title of the plot.
        Default is "Butte Creek Stream Transect".
    title_kwargs: dict, optional
        Additional keyword arguments for the title.
        Default is None.
    well_id_col: str, optional
        Column name for well IDs in obs_lut_df.
        Default is "site_code".
    well_annotations: list, optional
        List of dictionaries containing annotation parameters for well labels.
        If None, default annotations are created.
        Default is None.
    xlabel: str, optional
        Label for the x-axis.
        Default is "Distance (miles)".
    xlabel_kwargs: dict, optional
        Additional keyword arguments for the x-axis label.
        Default is None.
    xlim: tuple, optional
        Limits for the x-axis.
        If None, limits are determined from the data.
        Default is None.
    xticks_kwargs: dict, optional
        Additional keyword arguments for x-ticks.
        If None, default ticks are created.
        Default is None.
    ylim: tuple, optional
        Limits for the y-axis.
        If None, limits are determined from the data.
        Default is None.
    ylabel: str, optional
        Label for the y-axis.
        Default is "Elevation (ft amsl)".
    ylabel_kwargs: dict, optional
        Additional keyword arguments for the y-axis label.
        Default is None.


    Returns
    -------
    axes: matplotlib.axes.Axes
        The axes object containing the transect layout.
    fig: matplotlib.figure.Figure
        The figure object containing the transect layout.

    """

    borelog_id_col = kwargs.pop("borelog_id_col", "WELLINFOID")
    borelog_width = kwargs.pop("borelog_width", 1000)
    legend_items = kwargs.pop("legend_items", None)
    legend_kwargs = kwargs.pop("legend_kwargs", None)
    obs_col = kwargs.pop("obs_col", "gwe")
    obs_kwargs = kwargs.pop("obs_kwargs", None)
    proj_d_col = kwargs.pop("proj_d_col", "proj")
    sim_col = kwargs.pop("sim_col", "wte_ft")
    sim_kwargs = kwargs.pop("sim_kwargs", None)
    streambed_bottom_col = kwargs.pop("streambed_bottom_col", "botr")
    streambed_kwargs = kwargs.pop("streambed_kwargs", None)
    subplot_kwargs = kwargs.pop("subplot_kwargs", None)
    title = kwargs.pop("title", "Butte Creek Stream Transect")
    title_kwargs = kwargs.pop("title_kwargs", None)
    well_id_col = kwargs.pop("well_id_col", "site_code")
    well_annotations = kwargs.pop("well_annotations", None)
    xlabel = kwargs.pop("xlabel", "Distance (miles)")
    xlabel_kwargs = kwargs.pop("xlabel_kwargs", None)
    xlim = kwargs.pop("xlim", None)
    xticks_kwargs = kwargs.pop("xticks_kwargs", None)
    ylim = kwargs.pop("ylim", None)
    ylabel = kwargs.pop("ylabel", "Elevation (ft amsl)")
    ylabel_kwargs = kwargs.pop("ylabel_kwargs", None)


    if sim_kwargs is None:
        sim_kwargs = {
            "color": "navy",
            "label": "Simulated GWE"
        }

    if legend_kwargs is None:
        legend_kwargs = {
            "loc": "lower center",
            "bbox_to_anchor": (0.5, 0),
            "ncols": 5
        }

    if streambed_kwargs is None:
        streambed_kwargs = {
            "color": "#555555",
            "linestyle": "dashed",
            "label": "Streambed Bottom"
        }

    if obs_kwargs is None:
        obs_kwargs = {
            "label": 'Observed GWE',
            "markerfacecolor": "#4B7164",
            "marker": 'o',
            "color": "w"
        }

    if legend_items is None:
        legend_items = [
            Line2D(
                [0],
                [0],
                **obs_kwargs
            ),
            Line2D(
                [0],
                [0],
                **sim_kwargs
            ),
            Line2D(
                [0],
                [0],
                **streambed_kwargs
            )
        ]

    if subplot_kwargs is None:
        subplot_kwargs = {"nrows": 2,
                          "height_ratios": [0.9, 0.1],
                          "figsize": (12, 8)}

    if title_kwargs is None:
        title_kwargs = {"pad": 150}

    # Let's set default xlabel kwargs
    if xlabel_kwargs is None:
        xlabel_kwargs = {
        }

    # Let's make default xlim if none provided
    if xlim is None:
        xmin = stream_nodes_df["proj"].min()
        xmax = stream_nodes_df["proj"].max()
        xlim = (xmin, xmax)

    if ylabel_kwargs is None:
        ylabel_kwargs = {}

    if ylim is None:
        ymin = sim_df[sim_col].min()
        ymin = np.minimum(
            ymin,
            obs_df[obs_col].min()
        )
        ymin = np.minimum(
            ymin,
            stream_nodes_df[streambed_bottom_col].min()
        )
        ymax = sim_df[sim_col].max()
        ymax = np.maximum(
            ymax,
            obs_df[obs_col].max()
        )
        ymax = np.maximum(
            ymax,
            stream_nodes_df[streambed_bottom_col].max()
        )
        ylim = (ymin-10, ymax+10)

    # Make default well label annotations if none provided
    if well_annotations is None:
        well_annotations = []
        well_annotations_y = ylim[0] + 1.1 * (ylim[1] - ylim[0])
        annotation_kwargs = {
            "rotation": "vertical",
            "annotation_clip": False
        }
        for row in obs_lut_df.itertuples():
            proj_d = getattr(row, proj_d_col)
            annotation_kwargs_local = annotation_kwargs.copy()
            annotation_kwargs_local["text"] = getattr(row, well_id_col)
            annotation_kwargs_local["xy"] = (proj_d, ylim[1])
            annotation_kwargs_local["xytext"] = (proj_d, well_annotations_y)
            well_annotations.append(annotation_kwargs_local)

    # Let's make the default xticks kwargs
    if xticks_kwargs is None:
        step_miles = 5 * np.power(
            10, np.modf
                    (np.log10(
                    (
                            xlim[1] - xlim[0]
                    ) / 1609.34
                )
                )[1] - 1
        )
        ticks = np.arange(
            xlim[0],
            xlim[1],
            step_miles * 1609.34)
        labels = [f"{int(tick / 1609.34):,.0f}" for tick in ticks]
        xticks_kwargs = {
            "ticks": xlim[1] - ticks,
            "labels": labels
        }

    fig, axes = plt.subplots(**subplot_kwargs)

    # Well labels
    for annotation in well_annotations:
        axes[0].annotate(**annotation)

    # Lithology logs
    for row in lith_lut_df.itertuples():
        borelog_id = getattr(row, borelog_id_col)
        lith_df_local = lith_df.loc[
            lith_df[borelog_id_col] == borelog_id
            ].reset_index(drop=True)
        proj_d_borelog = getattr(row, proj_d_col)
        axes[0], legend_items = add_bore_logs(
            axes[0],
            lith_df_local,
            proj_d_borelog,
            patch_legend_items=legend_items,
            w=borelog_width
        )

    # Let's add streambed bottom elevations

    axes[0].plot(
        stream_nodes_df[proj_d_col],
        stream_nodes_df[streambed_bottom_col],
        **streambed_kwargs)

    # Set xlabel
    axes[0].set_xlabel(xlabel, **xlabel_kwargs)

    # Set ylabel
    axes[0].set_ylabel(ylabel, **ylabel_kwargs)

    # Set limits
    axes[0].set_xlim(xlim)
    axes[0].set_ylim(ylim)

    # Invert x-axis
    axes[0].xaxis.set_inverted(True)

    # Set x-ticks
    axes[0].set_xticks(**xticks_kwargs)

    # Set title
    axes[0].set_title(title, **title_kwargs)

    # Let's remove axis of the lower subplot
    axes[1].set_axis_off()

    # Let's add legend
    axes[1].legend(
        handles=legend_items,
        **legend_kwargs
    )

    plt.tight_layout()

    return fig, axes