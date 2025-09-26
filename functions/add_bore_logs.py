
import matplotlib.patches as patch

def add_bore_logs(
        ax,
        df,
        xdist,
        w=0.05,
        lith_top_col="LITH_TOP_DEPTH_ft",
        lith_bot_col="LITH_BOT_DEPTH_ft",
        lith_col="LITH_MAJOR1",
        styles=None,
        patch_legend_items = None,
        gse_col="GROUND_SURFACE_ELEVATION_ft"
):
    """
    Add AEM lithologies to a matplotlib axis object.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axis object to which the bore logs will be added.

    df : pandas.DataFrame
        Dataframe containing lithology information. Each row represents a horizon.
         Must contain the following columns:
        - 'GROUND_SURFACE_ELEVATION_ft' (float): Ground surface elevation in feet.
        - lith_top_col (str): Column name for the top depth of the lithology in feet.
        - lith_bot_col (str): Column name for the bottom depth of the lithology in feet.
        - lith_col (str): Column name for the lithology type.

    xdist : float
        X-coordinate where the bore log will be placed on the axis.

    w : float, optional
        Width of the bore log in x-axis units. Default is 0.05.

    lith_top_col : str, optional
        Column name for the top depth of the lithology in feet. Default is 'LITH_TOP_DEPTH_ft'.

    lith_bot_col : str, optional
        Column name for the bottom depth of the lithology in feet. Default is 'LITH_BOT_DEPTH_ft'.

    lith_col : str, optional
        Column name for the lithology type. Default is 'LITH_MAJOR1'.

    styles : dict, optional
        Dictionary defining the styles for different lithologies. Keys are lithology names, and values are
        dictionaries with 'facecolor' and 'hatch' keys. If None, default styles are used.

    patch_legend_items: list, optional
        List of matplotlib.patches.Rectangle objects for legend creation. If None, a new list will be created.
        Default is None.
    gse_col : str, optional
        Column name for the ground surface elevation in feet. Default is 'GROUND_SURFACE_ELEVATION_ft'.
        Default is 'GROUND_SURFACE_ELEVATION_ft'.

    Returns
    -------

    ax : matplotlib.axes.Axes
        The modified matplotlib axis object with the bore logs added.
    patch_legend_items : list
        List of matplotlib.patches.Rectangle objects for legend creation.

"""

    if styles is None:
        styles = {
            "clay": {"facecolor": "blueviolet", "hatch": "---///"},
            "silt": {"facecolor": "blue", "hatch": "xx"},
            "sand": {"facecolor": "yellow", "hatch": "..."},
            "gravel": {"facecolor": "orange", "hatch": "oo"},
            "cobbles": {"facecolor": "orange", "hatch": "oo"},
            "boulders": {"facecolor": "orange", "hatch": "oo"},
            "rock - undifferentiated": {"facecolor": "grey", "hatch": None},
            "rock - sedimentary": {"facecolor": "grey", "hatch": None},
            "rock - intrusive": {"facecolor": "grey", "hatch": None},
            "rock - metamorphic": {"facecolor": "grey", "hatch": None},
            "rock - volcanic": {"facecolor": "grey", "hatch": None},
            "soil": {"facecolor": "brown", "hatch": "/o"},
            "unknown": {"facecolor": "white", "hatch": "***"}
        }

    for i in range(len(df)):
        top_elevation = df.loc[i, gse_col]
        top_depth = df.loc[i, lith_top_col]
        height = top_depth - df.loc[i, lith_bot_col]
        color = styles[df.loc[i, lith_col]]["facecolor"]
        hatch = styles[df.loc[i, lith_col]]["hatch"]
        p = patch.Rectangle((xdist, top_elevation-top_depth), width=w, height=height, facecolor=color, edgecolor=None, hatch=hatch)
        ax.add_patch(p)

    if patch_legend_items is None:
        patch_legend_items = []
    for lith in ["clay", "silt", "sand", "soil", "unknown"]:
        if lith not in [patch_local.get_label().lower() for patch_local in patch_legend_items]:
            patch_legend_items.append(
                patch.Rectangle(
                    (0,0),
                    width=0,
                    height=0,
                    facecolor=styles[lith]["facecolor"],
                    edgecolor=None,
                    hatch=styles[lith]["hatch"],
                    label=lith.capitalize())
            )

    if any(x in df[lith_col].values for x in ["gravel", "cobbles", "boulders"]):
        if "gravel" not in [patch_local.get_label().lower() for patch_local in patch_legend_items]:
            patch_legend_items.append(
                patch.Rectangle(
                    (0,0),
                    width=0,
                    height=0,
                    facecolor=styles["gravel"]["facecolor"],
                    edgecolor=None,
                    hatch=styles["gravel"]["hatch"],
                    label="Gravel")
            )

    if any(x in df[lith_col].values for x in [
        "rock - undifferentiated",
        "rock - sedimentary",
        "rock - intrusive",
        "rock - metamorphic",
        "rock - volcanic"]
           ):
        if "rock" not in [patch_local.get_label().lower() for patch_local in patch_legend_items]:
            patch_legend_items.append(
                patch.Rectangle(
                    (0,0),
                    width=0,
                    height=0,
                    facecolor=styles["rock - undifferentiated"]["facecolor"],
                    edgecolor=None,
                    hatch=styles["rock - undifferentiated"]["hatch"],
                    label="Rock")
            )


    return ax, patch_legend_items
