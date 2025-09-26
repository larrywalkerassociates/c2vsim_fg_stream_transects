""" Function to make an animation of a transect figure containing borehole lithology logs, streambed bottom elevations,
 and simulated and observed groundwater elevations."""

import functools

from matplotlib import animation
import matplotlib.pyplot as plt
import pandas as pd

from functions.make_transect_layout import make_transect_layout


def make_transect_animation(
        stream_nodes_df,
        sim_df,
        obs_df,
        obs_lut_df,
        lith_lut_df,
        lith_df,
        **kwargs
)-> str:

    """
    Makes an animation of a transect figure containing borehole lithology logs, streambed bottom elevations,
    and simulated and observed groundwater elevations.

    Parameters
    ----------
    stream_nodes_df : pd.DataFrame
        Dataframe containing stream nodes information.
        Must contain columns with the following names (or names specified in keyword arguments):
        - 'igw': groundwater node index
        - 'proj': projection distance along the transect
    sim_df : pd.DataFrame
        Dataframe containing simulated groundwater elevations.
        Must contain columns with the following names (or names specified in keyword arguments):
        - 'date': date of the simulation
        - 'igw': groundwater node index
        - 'wte_ft': simulated groundwater elevation in feet
    obs_df : pd.DataFrame
        Dataframe containing observed groundwater elevations.
        Must contain columns with the following names (or names specified in keyword arguments):
        - 'date': date of the observation
        - 'site_code': well identifier
        - 'gwe': observed groundwater elevation in feet
    obs_lut_df : pd.DataFrame
        Dataframe containing lookup table for observed wells.
        Must contain columns with the following names (or names specified in keyword arguments):
        - 'site_code': well identifier
        - 'proj': projection distance along the transect
    lith_lut_df : pd.DataFrame
        Dataframe containing lookup table for lithology logs.
        Must contain columns with the following names (or names specified in keyword arguments):
        - 'WELLINFOID': log identifier
        - 'proj': projection distance along the transect
    lith_df : pd.DataFrame
        Dataframe containing lithology log information.
        Must contain columns with the following names (or names specified in keyword arguments):
        - 'WELLINFOID': log identifier
        - 'LITH_TOP_DEPTH_ft': top depth of the lithology unit
        - 'LITH_BOT_DEPTH_ft': bottom depth of the lithology unit
        - 'LITH_MAJOR1': lithology description for the lithology unit
        - 'GROUND_SURFACE_ELEVATION_ft': ground surface elevation at the log location

    Keyword Args
    ----------
    date_col: str, optional
        Column name for date in sim_df and obs_df.
        Default is 'date'.
    obs_col: str, optional
        Column name for observed groundwater elevation in obs_df.
        Default is 'gwe'.
    obs_kwargs: dict, optional
        Keyword arguments for customizing the appearance of observed data points.
        Default is None, which uses {'color': '#4B7164'}.
    proj_d_col: str, optional
        Column name for projection distance along the transect in stream_nodes_df, obs_lut_df
        and lith_lut_df.
        Default is 'proj'.
    sim_col: str, optional
        Column name for simulated groundwater elevation in sim_df.
        Default is 'wte_ft'.
    sim_kwargs: dict, optional
        Keyword arguments for customizing the appearance of simulated data line.
        Default is None, which uses {'color': 'navy', 'label': 'Simulated GWE'}.
    title: str, optional
        Title of the first frame of the animation.
        Default is 'Butte Creek Stream Transect'.
    title_kwargs: dict, optional
        Keyword arguments for customizing the appearance of the title.
        Default is None, which uses {'pad': 150}.
    well_id_col: str, optional
        Column name for well identifier in obs_df and obs_lut_df.
        Default is 'site_code'.
    node_key_col: str, optional
        Column name for groundwater node index in stream_nodes_df and sim_df.
        Default is 'igw'.
    interval: int, optional
        Interval between frames in milliseconds.
        Default is 500 milliseconds.


    Returns
    -------
    html_video: str
        HTML5 video string of the animation.

    """

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
        obs_kwargs = {"color": "#4B7164"}
        kwargs["obs_kwargs"] = obs_kwargs

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

    return html_video