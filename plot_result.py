import pandas as pd
import matplotlib.pyplot as plt

# Load datasets
df_LP = pd.read_csv("results_LP.csv")
df_QP = pd.read_csv("results_QP.csv")
df_MATLAB = pd.read_csv("results_MATLAB.csv")

datasets = [
    ("LP", df_LP),
    ("QP", df_QP),
    ("MATLAB MPC", df_MATLAB)
]

# Metrics
for name, df in datasets:
    mse = ((df['x1'] - df['x_r']) ** 2).mean()
    u_total = df['u'].abs().sum()

    print(f"{name} MSE (x1 vs setpoint): {mse:.4f}")
    print(f"{name} Total Control Effort (sum of |u|): {u_total:.4f}")
    print()

# Plot layout: 2 rows x 3 columns
fig, axes = plt.subplots(
    2,
    3,
    figsize=(18, 8),
    sharex=True
)

for col, (name, df) in enumerate(datasets):

    ax1 = axes[0, col]
    ax2 = axes[1, col]

    # States
    ax1.plot(
        df['k'],
        df['x0'],
        label='State 0 (x0)',
        color='orange',
        linestyle='--'
    )

    ax1.plot(
        df['k'],
        df['x1'],
        label='State 1 (x1)',
        color='blue',
        linewidth=2
    )

    ax1.plot(
        df['k'],
        df['x_r'],
        color='red',
        linestyle=':',
        label='Setpoint'
    )

    ax1.set_title(f'{name} — System States')
    ax1.set_ylabel('State Values')
    ax1.grid(True)
    ax1.legend()

    # Control effort
    ax2.step(
        df['k'],
        df['u'],
        label='Control Input (u)',
        color='green'
    )

    ax2.axhline(
        y=25.0,
        color='red',
        linestyle=':',
        label='Max u'
    )

    ax2.axhline(
        y=-25.0,
        color='purple',
        linestyle=':',
        label='Min u'
    )

    ax2.set_title(f'{name} — Control Effort')
    ax2.set_xlabel('Time Step (k)')
    ax2.set_ylabel('Control Effort')
    ax2.grid(True)
    ax2.legend()

plt.tight_layout()
plt.show()