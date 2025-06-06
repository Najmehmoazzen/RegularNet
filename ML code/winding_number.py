# %% [markdown]
# ## Regular Network


# %%
import numpy as np
import matplotlib.pyplot as plt

# ---------- تابع برای محاسبه winding number ----------
def compute_winding_number(theta_final):
    dtheta = np.angle(np.exp(1j * (theta_final - np.roll(theta_final, -1))))
    q = np.round(np.sum(dtheta) / (2 * np.pi))
    return int(q)

# ---------- اجرای چند شبیه‌سازی ----------

def run_trials(num_trials=50, output_file="winding_numbers.txt"):
    q_list = []
    with open(output_file, "w") as f:
        f.write("Sample\tWindingNumber\n")  # Header
        for trial in range(num_trials):
            filename = f"./kuramoto_cpp/Save/Last_Phase/Last_Phase_S_{trial}.txt"
            theta_final = np.loadtxt(filename)
            q = compute_winding_number(theta_final)
            q_list.append(q)
            f.write(f"{trial}\t{q}\n")
            #if q == -2:
            #    print(f"Sample {trial} has a winding number of {q}")
    return q_list

# %%
# ---------- اجرای شبیه‌سازی ----------
q_results = run_trials(num_trials=1700)

# ---------- نمایش نتایج ----------
unique_qs, counts = np.unique(q_results, return_counts=True)

plt.bar(unique_qs, counts, width=0.6)
plt.xlabel("Winding Number q")
plt.ylabel("Frequency of occurrence")
plt.title("Estimated basin sizes of different winding numbers")
plt.grid(True)
plt.savefig("dist.jpg", dpi=100)



