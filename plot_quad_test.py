"""Plot results of tests/mpc_quad_hover_test.cpp.

Usage (from the folder containing the CSVs):
    python plot_quad_mpc_test.py
    python plot_quad_mpc_test.py --dir build --save
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

G = 9.81
TILT_MAX = 0.35  # [rad]
T_MIN = 0.3 * G  # [m/s^2]
T_MAX = 1.8 * G  # [m/s^2]

SCENARIOS = {
    'Hover hold': 'results_quad_hover.csv',
    'Altitude step': 'results_quad_alt.csv',
    'Lateral step': 'results_quad_lat.csv',
    'Mid-run lateral step': 'results_quad_lat_mid.csv',
}


def plot_scenario(df, title):
    """Plot states, inputs, solve time and iterations of one scenario."""
    fig, axes = plt.subplots(6, 1, figsize=(10, 14), sharex=True)
    t = df['t']

    ax = axes[0]
    for name in ('px', 'py', 'pz'):
        line, = ax.plot(t, df[name], label=name)
        ref = df[f'{name}_ref']
        if ref.abs().max() > 0.0:
            ax.step(t, ref, where='post', color=line.get_color(),
                    linestyle=':', label=f'{name} ref')
    ax.set_ylabel('Position [m]')

    ax = axes[1]
    for name in ('vx', 'vy', 'vz'):
        ax.plot(t, df[name], label=name)
    ax.set_ylabel('Velocity [m/s]')

    ax = axes[2]
    ax.plot(t, df['phi'], label='phi')
    ax.plot(t, df['theta'], label='theta')
    ax.step(t, df['phi_r'], where='post', linestyle='--', label='phi_r')
    ax.step(t, df['theta_r'], where='post', linestyle='--', label='theta_r')
    ax.axhline(TILT_MAX, color='red', linestyle=':', label='limit')
    ax.axhline(-TILT_MAX, color='red', linestyle=':')
    ax.set_ylabel('Attitude [rad]')

    ax = axes[3]
    ax.step(t, df['T'], where='post', label='T')
    ax.axhline(G, color='gray', linestyle='--', label='hover')
    ax.axhline(T_MIN, color='red', linestyle=':', label='limit')
    ax.axhline(T_MAX, color='red', linestyle=':')
    ax.set_ylabel('Thrust [m/s^2]')

    ax = axes[4]
    ax.bar(t, df['solve_ms'], width=0.08, label='solve time')
    ax.axhline(df['solve_ms'].mean(), color='black', linestyle='--',
               label=f"avg {df['solve_ms'].mean():.2f} ms")
    ax.set_ylabel('Solve [ms]')

    ax = axes[5]
    ax.bar(t, df['iters'], width=0.08, color='tab:orange', label='iterations')
    ax.axhline(df['iters'].mean(), color='black', linestyle='--',
               label=f"avg {df['iters'].mean():.1f}")
    ax.set_ylabel('Iterations')
    ax.set_xlabel('Time [s]')

    for ax in axes:
        ax.grid(True)
        ax.legend(loc='upper right', fontsize='small')

    fig.suptitle(title)
    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dir', default='.', help='folder with the CSVs')
    parser.add_argument('--save', action='store_true',
                        help='save PNGs instead of showing the plots')
    args = parser.parse_args()
    data_dir = Path(args.dir)

    print(f"{'Scenario':<22} {'avg [ms]':>10} {'worst [ms]':>11} "
          f"{'avg iters':>10} {'worst iters':>12} {'ms/iter':>8}")
    for title, filename in SCENARIOS.items():
        path = data_dir / filename
        if not path.exists():
            print(f'{title:<22} missing: {path}')
            continue

        df = pd.read_csv(path)
        worst = df['iters'].idxmax()
        ms_per_iter = df['solve_ms'][worst] / max(df['iters'][worst], 1)
        print(f"{title:<22} {df['solve_ms'].mean():>10.3f} "
              f"{df['solve_ms'].max():>11.3f} {df['iters'].mean():>10.2f} "
              f"{df['iters'].max():>12d} {ms_per_iter:>8.3f}")

        fig = plot_scenario(df, title)
        if args.save:
            fig.savefig(data_dir / path.with_suffix('.png').name, dpi=120)
            plt.close(fig)

    if not args.save:
        plt.show()


if __name__ == '__main__':
    main()