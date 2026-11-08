"""Plot the trajectory scenario of tests/linear_drone.cpp.

README layout: xy path on top, roll and pitch below, solve time per step
with the real-time target at the bottom.

Usage (from the folder containing the CSV):
    python plot_quad_trajectory.py
    python plot_quad_trajectory.py --save
    python plot_quad_trajectory.py --file build/results_quad_trajectory.csv
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TILT_MAX_DEG = np.degrees(0.35)
SOLVE_TARGET_MS = 10.0

# Gates of the trajectory scenario, must match linear_drone.cpp
GATES = np.array([[3.0, 2.0], [6.0, -2.0], [9.0, 2.0],
                  [12.0, -2.0], [15.0, 0.0]])


def plot_xy(ax, df):
    ax.plot(df['px_ref'], df['py_ref'], linestyle=':', marker='.',
            markersize=3, color='k', alpha=0.6, zorder=2,
            label='Reference trajectory (every 0.1 s)')
    ax.plot(df['px'], df['py'], linewidth=2, zorder=3, label='Drone path')
    ax.scatter(GATES[:, 0], GATES[:, 1], marker='x', color='tab:red', s=60,
               zorder=4, label='Gates')
    ax.scatter(df['px'].iloc[0], df['py'].iloc[0], color='g', zorder=4,
               label='Start')
    ax.set_xlabel('x position [m]')
    ax.set_ylabel('y position [m]')
    ax.set_title('Flight path (top view)')
    ax.set_aspect('equal')


def plot_angle(ax, t, state, ref, name, symbol, direction):
    ax.plot(t, np.degrees(state), label=f'{name} angle {symbol} (state)')
    ax.step(t, np.degrees(ref), where='post', linestyle='--',
            label=f'{name} command {symbol}_r (MPC output)')
    for sign in (1, -1):
        ax.axhline(sign * TILT_MAX_DEG, color='r', linestyle=':',
                   label='Tilt limit ±20°' if sign == 1 else None)
    ax.set_title(f'{name}: tilts thrust along {direction}')
    ax.set_ylabel(f'{name} {symbol} [deg]')
    ax.set_xlabel('Time [s]')


def plot_solve_time(ax, df):
    worst = df['solve_ms'].max()
    ax.step(df['t'], df['solve_ms'], where='post',
            label=f'QP solve time (worst {worst:.1f} ms)')
    ax.axhline(SOLVE_TARGET_MS, color='r', linestyle='--',
               label=f'Real-time target {SOLVE_TARGET_MS:.0f} ms')
    ax.set_title('Active-set QP solver, one solve per 0.1 s control step')
    ax.set_ylabel('Solve time [ms]')
    ax.set_xlabel('Time [s]')
    ax.set_ylim(bottom=0.0)


def plot_scenario(df, title):
    fig = plt.figure(figsize=(13, 10), layout='constrained')
    grid = fig.add_gridspec(3, 2, height_ratios=(1.4, 1, 0.8))
    ax_xy = fig.add_subplot(grid[0, :])
    ax_roll = fig.add_subplot(grid[1, 0])
    ax_pitch = fig.add_subplot(grid[1, 1], sharex=ax_roll, sharey=ax_roll)
    ax_time = fig.add_subplot(grid[2, :])
    t = df['t']

    plot_xy(ax_xy, df)
    plot_angle(ax_roll, t, df['phi'], df['phi_r'], 'Roll', 'φ', '-y')
    plot_angle(ax_pitch, t, df['theta'], df['theta_r'], 'Pitch', 'θ', '+x')
    plot_solve_time(ax_time, df)

    for ax in (ax_xy, ax_roll, ax_pitch, ax_time):
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='upper right')
    fig.suptitle(title, fontsize=14)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--file', type=Path,
                        default=Path('results_quad_trajectory.csv'),
                        help='trajectory CSV written by linear_drone')
    parser.add_argument('--save', action='store_true',
                        help='save a PNG next to the CSV instead of showing')
    args = parser.parse_args()

    df = pd.read_csv(args.file)
    fig = plot_scenario(df, 'Linear MPC: trajectory tracking through gates')
    if args.save:
        out = args.file.with_suffix('.png')
        fig.savefig(out, dpi=150)
        print(f'saved {out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()