"""Plot the results of tests/pid_test.cpp.

Compares the velocity form and positional PID controllers on the same
2-state system and reference as tests/mpc_test.cpp. If results_QP.csv from
mpc_test is in the same folder, the QP-MPC result is drawn for comparison.

Usage (from the repo root, after ./build/pid_test):
    python plot_pid_test.py
    python plot_pid_test.py --save
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

U_LIMIT = 25.0

RUNS = {
    'PID velocity form': 'results_PID_velocity.csv',
    'PID positional form': 'results_PID_positional.csv',
    'MPC (QP cost)': 'results_QP.csv',
}


def load_runs(data_dir):
    """Return {label: DataFrame} for every result file that exists."""
    runs = {}
    for label, name in RUNS.items():
        path = data_dir / name
        if path.exists():
            runs[label] = pd.read_csv(path)
        elif 'PID' in label:
            raise FileNotFoundError(f'{path} not found, run pid_test first')
    return runs


def plot_runs(runs):
    fig, (ax_y, ax_u) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                     layout='constrained')
    first = next(iter(runs.values()))
    ax_y.step(first['k'], first['x_r'], where='post', color='k',
              linestyle=':', label='Reference')

    for label, df in runs.items():
        ax_y.plot(df['k'], df['x1'], label=f'{label}: output x1')
        ax_u.step(df['k'], df['u'], where='post', label=f'{label}: input u')

    for sign in (1, -1):
        ax_u.axhline(sign * U_LIMIT, color='r', linestyle=':',
                     label=f'Input limit ±{U_LIMIT:.0f}' if sign == 1
                     else None)

    ax_y.set_ylabel('Output x1')
    ax_y.set_title('Controlled output')
    ax_u.set_ylabel('Control input u')
    ax_u.set_xlabel('Time step k')
    ax_u.set_title('Control input')
    for ax in (ax_y, ax_u):
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle('PID controllers on the 2-state test system')
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dir', type=Path, default=Path('.'),
                        help='folder with the CSV files')
    parser.add_argument('--save', action='store_true',
                        help='save pid_results.png instead of showing')
    args = parser.parse_args()

    runs = load_runs(args.dir)
    for label, df in runs.items():
        mse = ((df['x1'] - df['x_r']) ** 2).mean()
        effort = df['u'].abs().sum()
        print(f'{label:22s} MSE {mse:8.4f}   sum|u| {effort:9.2f}')

    fig = plot_runs(runs)
    if args.save:
        out = args.dir / 'pid_results.png'
        fig.savefig(out, dpi=150)
        print(f'saved {out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()