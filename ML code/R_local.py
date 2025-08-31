import numpy as np
import os

def local_R(sample, data):
    N = len(data)         # Total number of nodes
    k = 2                 # 2 neighbors on each side
    R = np.zeros(N)       # Array to store local order parameter for each node
    
    for i in range(N):
        # Neighborhood indices with periodic boundary conditions (ring)
        neighborhood = [(i + j) % N for j in range(-k, k+1)]
        local_phases = data[neighborhood]
        local_order = np.abs(np.mean(np.exp(1j * local_phases)))
        R[i] = local_order

    # Save only R (one column)
    np.savetxt(
        f"./local_R/local_R_{sample}.txt",
        R,
        fmt="%.4f",
        comments=""
    )
    return R


# === Process all samples ===
num_samples = 40000
os.makedirs("./local_R", exist_ok=True)  # Ensure output folder exists

for sample in range(num_samples):
    data = np.loadtxt(f"./kuramoto_cpp/Save/I=InitialPhases/S_{sample}.txt")  # Load phases
    local_R(sample, data)  # Compute and save local order parameter
    
    # Print progress every 1000 samples
    if sample % 1000 == 0:
        print(f"Processed sample {sample}/{num_samples}")
