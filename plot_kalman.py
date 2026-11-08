"""Plot Kalman filter results from kalman_results.csv."""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSV_PATH = sys.argv[1] if len(sys.argv) > 1 else "kalman_results.csv"

df = pd.read_csv(CSV_PATH)

pos_err = df["x_est_0"] - df["x_true_0"]
vel_err = df["x_est_1"] - df["x_true_1"]
sigma_pos = np.sqrt(df["P00"])
sigma_vel = np.sqrt(df["P11"])

fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

# Position: measurement is noisy, estimate should sit close to truth.
ax = axes[0]
ax.plot(df["k"], df["z"], ".", color="0.7", label="measurement")
ax.plot(df["k"], df["x_true_0"], "k-", label="true")
ax.plot(df["k"], df["x_est_0"], "r--", label="estimate")
ax.set_ylabel("position")
ax.legend(loc="upper left")
ax.grid(alpha=0.3)

# Velocity is unmeasured, so this panel is the real test.
ax = axes[1]
ax.plot(df["k"], df["x_true_1"], "k-", label="true")
ax.plot(df["k"], df["x_est_1"], "r--", label="estimate")
ax.set_ylabel("velocity")
ax.legend(loc="upper left")
ax.grid(alpha=0.3)

# Errors against the filter's own 2-sigma band: roughly 95% should fall inside.
ax = axes[2]
ax.plot(df["k"], pos_err, label="position error")
ax.plot(df["k"], vel_err, label="velocity error")
ax.fill_between(df["k"], -2 * sigma_pos, 2 * sigma_pos, alpha=0.15,
                color="C0", label=r"$\pm 2\sigma$ position")
ax.fill_between(df["k"], -2 * sigma_vel, 2 * sigma_vel, alpha=0.15,
                color="C1", label=r"$\pm 2\sigma$ velocity")
ax.axhline(0.0, color="k", lw=0.5)
ax.set_ylabel("error")
ax.legend(loc="upper left", fontsize=8)
ax.grid(alpha=0.3)

# Covariance should converge to a constant within 10-20 steps.
ax = axes[3]
ax.semilogy(df["k"], df["P00"], label="P00")
ax.semilogy(df["k"], df["P11"], label="P11")
ax.set_ylabel("covariance")
ax.set_xlabel("step k")
ax.legend(loc="upper right")
ax.grid(alpha=0.3, which="both")

fig.tight_layout()
fig.savefig("kalman_results.png", dpi=150)

tail = slice(-20, None)
inside_pos = (pos_err.abs() < 2 * sigma_pos).mean()
inside_vel = (vel_err.abs() < 2 * sigma_vel).mean()

print(f"steady-state P00: {df['P00'].iloc[tail].mean():.6f} "
      f"(std {df['P00'].iloc[tail].std():.2e})")
print(f"steady-state P11: {df['P11'].iloc[tail].mean():.6f} "
      f"(std {df['P11'].iloc[tail].std():.2e})")
print(f"position error inside 2-sigma: {inside_pos:.1%}")
print(f"velocity error inside 2-sigma: {inside_vel:.1%}")
print(f"position error mean: {pos_err.mean():+.4f}")
print(f"velocity error mean: {vel_err.mean():+.4f}")

plt.show()