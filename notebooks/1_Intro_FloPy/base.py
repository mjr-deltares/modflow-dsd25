import flopy
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.animation import FuncAnimation
import numpy as np

MOSAIC = [
    ["a", "a", "a", "c"],
    ["a", "a", "a", "c"],
    ["a", "a", "a", "c"],
    ["b", "b", "b", "."],
    ]
CROSS_SECTION_COLUMN = 25
CROSS_SECTION_ROW = 24

FPS = 4
ANI_EXT = ".mp4"
Writer = mpl.animation.writers["ffmpeg"]
writer = Writer(fps=FPS, metadata=dict(artist="dsd25"), bitrate=2056)


def plot_model_results(
        model, 
        cobj, 
        budobj, 
        totim=3652.5, 
        vmin=None, 
        vmax=None, 
        layer=0,
        path=None,
        dpi=300,
        ):
    if vmax is None:
        vmax = totim / 365.25
    if vmin is None:
        vmin = 1.0

    # get the simulated data for timestep
    c = cobj.get_data(totim=totim)
    spd = budobj.get_data(totim=totim, text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.get_specific_discharge(spd, model=model)

    b_y = model.modelgrid.ycellcenters[CROSS_SECTION_ROW, 0]
    c_x = model.modelgrid.xcellcenters[0, CROSS_SECTION_COLUMN]

    norm = None #colors.LogNorm(vmin=vmin, vmax=vmax)
    
    with flopy.plot.styles.USGSMap():
        fig, axd = plt.subplot_mosaic(
            MOSAIC,
            figsize=(10,8),
            layout="constrained",
            )
        
        fig.suptitle(f"Simulation time: {totim/365.25:5.2f} years - Layer {layer + 1}")
        
        ax = axd["a"]
        mv = flopy.plot.PlotMapView(model=model, ax=ax, layer=layer)
        pa = mv.plot_array(c, norm=norm, vmin=vmin, vmax=vmax)
        mv.plot_bc(package=model.river, color="blue", alpha=0.5)
        mv.plot_vector(qx, qy)
        ax.axhline(b_y, lw=2.0, color="red", ls=":")
        ax.axvline(c_x, lw=2.0, color="red", ls=":")
        ax.set_xticklabels([])
        ax.set_ylabel("y-position, m")

        axins = inset_axes(
            ax,
            width="25%",  # Width of the colorbar (as a percentage of parent axes width)
            height="2.5%", # Height of the colorbar (as a percentage of parent axes height)
            loc="lower left", # Location within the parent axes (e.g., 'upper right', 'lower left', 'center right')
            bbox_to_anchor=(0.04, 0.035, 1, 1), # Fine-tune position relative to parent's bbox
            bbox_transform=ax.transAxes,
            borderpad=0,
            )
        cbar = plt.colorbar(pa, cax=axins, orientation="horizontal")  
        cbar.ax.set_title("Groundwater age, years", fontdict={"fontsize": 8, "fontweight": "bold"})
      

        ax = axd["b"]
        xs = flopy.plot.PlotCrossSection(model=model, ax=ax, line={"row": CROSS_SECTION_ROW})
        pa = xs.plot_array(c, norm=norm, vmin=vmin, vmax=vmax)
        xs.plot_grid(lw=0.5, color="black")
        ax.axvline(c_x, lw=2.0, color="red", ls=":")
        ax.set_ylabel("Elevation, m")
        ax.set_xlabel("x-position, m")

        ax = axd["c"]
        xs = flopy.plot.PlotCrossSection(
            model=model, 
            ax=ax, 
            geographic_coords=False, 
            line={"column": CROSS_SECTION_COLUMN},
            )
        pa = xs.plot_array(c, norm=norm, vmin=vmin, vmax=vmax)
        pg = xs.plot_grid(lw=0.5, color="black")

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xc, yc = 0.5 * sum(xlim), 0.5 * sum(ylim)
        rotation_point = np.array([xc, yc])
        angle_deg = -90

        transform = (
            mtransforms.Affine2D().translate(*-rotation_point) +  # Translate to origin
            mtransforms.Affine2D().rotate_deg(angle_deg) +        # Rotate
            mtransforms.Affine2D().translate(*np.flip(rotation_point))      # Translate back
        )

        pa.set_transform(transform + ax.transData)
        pg.set_transform(transform + ax.transData)

        ax.set_xlim(ylim[::-1])
        ax.set_ylim(xlim)
        ax.axhline(b_y, lw=2.0, color="red", ls=":")
        ax.set_yticklabels([])
        ax.set_xlabel("Elevation, m")

        if path is not None:
            path.parent.mkdir(exist_ok=True, parents=True)
            fig.savefig(path, dpi=dpi, transparent=False)


def animate_model_results(
        model, 
        cobj, 
        budobj, 
        vmin, 
        vmax, 
        path,
        layer=0,
        ):

    path.parent.mkdir(exist_ok=True, parents=True)

    times = cobj.get_times()
    totim = times[0]
    ntimes = len(times)
    frames = np.arange(1, ntimes, 1)

    # get the simulated data for timestep
    c = cobj.get_data(totim=totim)
    spd = budobj.get_data(totim=totim, text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.get_specific_discharge(spd, model=model)

    b_y = model.modelgrid.ycellcenters[CROSS_SECTION_ROW, 0]
    c_x = model.modelgrid.xcellcenters[0, CROSS_SECTION_COLUMN]
    
    with flopy.plot.styles.USGSMap():
        fig, axd = plt.subplot_mosaic(
            MOSAIC,
            figsize=(10,8),
            layout="constrained",
            )
        
        fig.suptitle(f"Simulation time: {totim/365.25:3.0f} years - Layer {layer + 1}")
        
        ax = axd["a"]
        mv = flopy.plot.PlotMapView(model=model, ax=ax, layer=layer)
        pa_a = mv.plot_array(c, vmin=vmin, vmax=vmax)
        mv.plot_bc(package=model.river, color="blue", alpha=0.5)
        Q_a = mv.plot_vector(qx, qy)
        ax.axhline(b_y, lw=2.0, color="red", ls=":")
        ax.axvline(c_x, lw=2.0, color="red", ls=":")
        ax.set_xticklabels([])
        ax.set_ylabel("y-position, m")

        axins = inset_axes(
            ax,
            width="25%",  # Width of the colorbar (as a percentage of parent axes width)
            height="2.5%", # Height of the colorbar (as a percentage of parent axes height)
            loc="lower left", # Location within the parent axes (e.g., 'upper right', 'lower left', 'center right')
            bbox_to_anchor=(0.04, 0.035, 1, 1), # Fine-tune position relative to parent's bbox
            bbox_transform=ax.transAxes,
            borderpad=0,
            )
        cbar = plt.colorbar(pa_a, cax=axins, orientation="horizontal")  
        cbar.ax.set_title("Groundwater age, years", fontdict={"fontsize": 8, "fontweight": "bold"})
      

        ax = axd["b"]
        xs = flopy.plot.PlotCrossSection(model=model, ax=ax, line={"row": CROSS_SECTION_ROW})
        pa_b = xs.plot_array(c, vmin=vmin, vmax=vmax)
        xs.plot_grid(lw=0.5, color="black")
        ax.axvline(c_x, lw=2.0, color="red", ls=":")
        ax.set_ylabel("Elevation, m")
        ax.set_xlabel("x-position, m")

        ax = axd["c"]
        xs = flopy.plot.PlotCrossSection(
            model=model, 
            ax=ax, 
            geographic_coords=False, 
            line={"column": CROSS_SECTION_COLUMN},
            )
        pa_c = xs.plot_array(c, vmin=vmin, vmax=vmax)
        pg = xs.plot_grid(lw=0.5, color="black")

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xc, yc = 0.5 * sum(xlim), 0.5 * sum(ylim)
        rotation_point = np.array([xc, yc])
        angle_deg = -90

        transform = (
            mtransforms.Affine2D().translate(*-rotation_point) +  # Translate to origin
            mtransforms.Affine2D().rotate_deg(angle_deg) +        # Rotate
            mtransforms.Affine2D().translate(*np.flip(rotation_point))      # Translate back
        )

        pa_c.set_transform(transform + ax.transData)
        pg.set_transform(transform + ax.transData)

        ax.set_xlim(ylim[::-1])
        ax.set_ylim(xlim)
        ax.axhline(b_y, lw=2.0, color="red", ls=":")
        ax.set_yticklabels([])
        ax.set_xlabel("Elevation, m")

        def func(idx):

            totim = times[idx]
            fig.suptitle(f"Simulation time: {totim/365.25:5.2f} years - Layer {layer + 1}")

            c = cobj.get_data(totim=totim)
            c[c == 1e30] = np.nan
            c_b = c[:, CROSS_SECTION_ROW, :].flatten()
            c_b = c_b[~np.isnan(c_b)]
            c_c = c[:, :, CROSS_SECTION_COLUMN].flatten()
            c_c = c_c[~np.isnan(c_c)]
            spd = budobj.get_data(totim=totim, text="DATA-SPDIS")[0]
            qx, qy, _ = flopy.utils.get_specific_discharge(spd, model=model)
    
            pa_a.set_array(c[layer, :, :])
            pa_b.set_array(c_b)
            pa_c.set_array(c_c)

            Q_a.set_UVC(qx[layer, :, :], qy[layer, :, :])

            return pa_a, pa_b, pa_c, Q_a
        
        ani = FuncAnimation(fig, func, frames=frames, blit=False)
        plt.close()

    path = path.with_suffix(ANI_EXT)
    ani.save(path, writer=writer)
