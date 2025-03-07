"""
Plot errors for the raw and truncated series against JPL DE 441.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

from testing import JPL_DE_ROUTES
import tools.misc as misc
from vsop87a_ephemeris import VSOP87AEphemeris
from mpp02_ephemeris import MPP02Ephemeris

VSOP87A_RAW_JSON_PATH = R'./json/vsop87a_raw.json'
VSOP87A_TRUNCATED_JSON_PATH = lambda size: f'./json/vsop87a_truncated_{size}.json'
MPP02_LLR_RAW_JSON_PATH = R'./json/mpp02_llr_raw.json'
MPP02_LLR_TRUNCATED_JSON_PATH = lambda size: f'./json/mpp02_llr_truncated_{size}.json'
# MPP02_405_RAW_JSON_PATH = R'./json/mpp02_405_raw.json'
JPL_DE_EPHEMERIS_PATH = R'd:/resources/astro/de/de441.bsp'

PLOT_SAVE_DIRECTORY = R'./images'

# ecliptic to equatorial J2000.0
ROT_ECL_EQU = misc.rotation_matrix(0, 84381.448/3600*np.pi/180.0)

def _moving_average_shrinking(data, w):
    n = len(data)
    w = min(w, (n+1)//2)
    cs = np.cumsum(data)
    y = data.copy()
    for k in range(1, w):
        y[k] = cs[2*k] / (2*k+1)
        y[n-k] = (cs[n-1] - cs[n-2*k]) / (2*k-1)
    if n-w+1 > w:
        # middle section
        cs_left = cs[:n-2*w+1]
        cs_right = cs[2*w-1:]
        y[w:n-w+1] = (cs_right - cs_left) / (2*w-1)
    return y

def _moving_average_nonshrinking(data, w):
    n = len(data)
    w = min(w, n)

    def avg(k):
        return np.average(data[max(k-w+1, 0):min(k+w, n)])
    
    cs = np.cumsum(data)
    y = np.zeros(n)
    for k in range(w):
        # lost patience with indexing cs - just brute force it
        y[k] = avg(k)
        y[n-k-1] = avg(n-k-1)
    if n-w+1 > w:
        # middle section
        cs_left = cs[:n-2*w+1]
        cs_right = cs[2*w-1:]
        y[w:n-w+1] = (cs_right - cs_left) / (2*w-1)
    return y

def moving_average(data, w, shrinking=False):
    if shrinking:
        return _moving_average_shrinking(data, w)
    return _moving_average_nonshrinking(data, w)

def plot_errors(jpl_pv, my_pv, title):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21, 6))
    colors = ['red', 'navy', 'green', 'dimgray']
    line_widths = [8, 4, 3, 2]
    timescales = [100.0, 20.0, 4.0]
    for ax_index, ax in enumerate([ax1, ax2, ax3]):
        t_min = -timescales[ax_index]
        t_max = timescales[ax_index]
        for index, size in enumerate(['raw', 'large', 'medium', 'small']):
            color = colors[index]
            line_width = line_widths[index]
            errors = []
            num = 20000
            for t in np.linspace(t_min, t_max, num):
                p_ref, _ = jpl_pv(t)
                p, _ = my_pv(size, t)
                err = np.linalg.norm(p-p_ref) / np.linalg.norm(p_ref)
                errors.append((2000.0+t*100.0, err/misc.ARCSEC))
            errors = np.array(errors)

            errors_avg = moving_average(errors[:,1], 100)

            ax.scatter(errors[:,0], errors[:,1], s=1, color=color, alpha=0.2)
            ax.plot(errors[:,0], errors_avg, color=color, linewidth=line_width, label=size if size != 'raw' else 'no truncation')

        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: str(int(val) if val.is_integer() else val)))
        # ax.set_xscale('symlog')
        ax.set_xlabel('Time [Julian year]', fontsize=14)
        if ax_index == 0:
            ax.set_ylabel(R'Relative error [arcsec]', fontsize=14)
        # ax.set_title(title)
        ax.set_axisbelow(True)
        ax.grid()
        ax.legend()

    fig.subplots_adjust(wspace=0.15)
    fig.suptitle(title, fontsize=18)

    return fig

def save_image(fig: plt.Figure, file_path):
    fig.savefig(file_path, dpi=80, bbox_inches='tight')
    print(f'Wrote image to file {file_path}.')
    plt.close(fig)

@misc.time_it
def main():
    sizes = ['raw', 'small', 'medium', 'large']
    mpp02 = {}
    for size in sizes:
        json_path = MPP02_LLR_TRUNCATED_JSON_PATH(size) if size != 'raw' else MPP02_LLR_RAW_JSON_PATH
        mpp02[size] = MPP02Ephemeris(json_path)
    vsop87a = {}
    for size in sizes:
        json_path = VSOP87A_TRUNCATED_JSON_PATH(size) if size != 'raw' else VSOP87A_RAW_JSON_PATH
        vsop87a[size] = VSOP87AEphemeris(json_path)

    with misc.jplephem_pos_vel(JPL_DE_EPHEMERIS_PATH) as jpl_pos_vel:
        # Moon
        jpl_pv_moon = lambda t: jpl_pos_vel([399, 3, 301], t)
        pv = lambda size, t: mpp02[size].get_pos_vel(t)
        fig = plot_errors(jpl_pv_moon, pv, 'MOON - Relative error of ELP/MPP02 and its truncations w.r.t. JPL DE441')
        save_image(fig, file_path = f'{PLOT_SAVE_DIRECTORY}/error_moon.jpg')

        # Planets
        for body_name in vsop87a['raw'].bodies.keys():
            if body_name == 'EARTH':
                continue
            jpl_pv_moon = lambda t: jpl_pos_vel(JPL_DE_ROUTES[body_name], t)
            pv = lambda size, t: vsop87a[size].get_pos_vel(body_name, t)
            fig = plot_errors(jpl_pv_moon, pv, f'{body_name} - Relative error of VSOP87A and its truncations w.r.t. JPL DE441')
            save_image(fig, file_path = f'{PLOT_SAVE_DIRECTORY}/error_{body_name.lower()}.jpg')

if __name__ == '__main__':
    main()