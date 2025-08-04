import numpy as np 
import matplotlib.pyplot as plt

# === 参数设定 ===
P_mean_ppm = (785 + 2032) / 2
P_std_ppm = (2032 - 785) / 4
P_release_fraction = 0.55
density = 2.23e3  # kg/m³
width_m = 300e3   # m
k = 0.01 / 1e3    # m⁻¹
n_samples = 10000

# === C/P（gC/gP）===
# C/P_molar = (182, 831)
cp_mass_min = 70.6
cp_mass_max = 322.3

# === space===
T0_cm_array = np.arange(0.001, 1.7, 0.002)
T0_m_array = T0_cm_array / 100

x_near_m = np.linspace(0, 150e3, 300)
x_far_m = np.linspace(150e3, 300e3, 300)
dx_near = x_near_m[1] - x_near_m[0]
dx_far = x_far_m[1] - x_far_m[0]

# === save ===
volume_near_km3 = []
volume_far_km3 = []

release_mean_near = []
release_std_near = []
release_mean_far = []
release_std_far = []

OC_mean_total = []  # TOC
OC_std_total = []   # TOC_std

# === main ===
for T0 in T0_m_array:

    T_near = T0 * np.exp(-k * x_near_m)
    V_near_m3 = np.sum(T_near * dx_near * width_m)
    V_near_km3 = V_near_m3 / 1e9
    M_near_kg = V_near_m3 * density

    T_far = T0 * np.exp(-k * x_far_m)
    V_far_m3 = np.sum(T_far * dx_far * width_m)
    V_far_km3 = V_far_m3 / 1e9
    M_far_kg = V_far_m3 * density

    volume_near_km3.append(V_near_km3)
    volume_far_km3.append(V_far_km3)

    P_samples = np.random.normal(P_mean_ppm, P_std_ppm, n_samples) / 1e6
    P_samples[P_samples < 0] = 0


    P_release_near = M_near_kg * P_samples * P_release_fraction / 1e3
    P_release_far = M_far_kg * P_samples * P_release_fraction / 1e3
    P_release_total = P_release_near + P_release_far  # Lake

    release_mean_near.append(np.mean(P_release_near))
    release_std_near.append(np.std(P_release_near))
    release_mean_far.append(np.mean(P_release_far))
    release_std_far.append(np.std(P_release_far))

    
    cp_mass_ratio = np.random.uniform(cp_mass_min, cp_mass_max, n_samples)

    OC_g_total = P_release_total * 1e6 * cp_mass_ratio  # gC
    OC_ton_total = OC_g_total / 1e6                     # ton_C

    OC_mean_total.append(np.mean(OC_ton_total))
    OC_std_total.append(np.std(OC_ton_total))


volume_near_km3 = np.array(volume_near_km3)
volume_far_km3 = np.array(volume_far_km3)
release_mean_near = np.array(release_mean_near)
release_std_near = np.array(release_std_near)
release_mean_far = np.array(release_mean_far)
release_std_far = np.array(release_std_far)
OC_mean_total = np.array(OC_mean_total)
OC_std_total = np.array(OC_std_total)

# === Fig1 ===
plt.figure(figsize=(9, 6))
plt.plot(T0_cm_array, volume_near_km3, label='Near (0–150 km)', color='blue')
plt.plot(T0_cm_array, volume_far_km3, label='Far (150–300 km)', color='green')
plt.xlabel("Initial Ash Thickness (cm)")
plt.ylabel("Ash Volume (km³)")
plt.title("Ash Volume vs Initial Thickness")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("ash_volume_vs_thickness.jpg", dpi=300)

# === Fig2 ===
plt.figure(figsize=(9, 6))
sigma = 1
plt.fill_between(T0_cm_array,
                 release_mean_near - sigma * release_std_near,
                 release_mean_near + sigma * release_std_near,
                 color='skyblue', alpha=0.5, label='Near ±1σ')
plt.plot(T0_cm_array, release_mean_near, color='blue', label='Near Mean')

plt.fill_between(T0_cm_array,
                 release_mean_far - sigma * release_std_far,
                 release_mean_far + sigma * release_std_far,
                 color='lightgreen', alpha=0.5, label='Far ±1σ')
plt.plot(T0_cm_array, release_mean_far, color='green', label='Far Mean')

plt.xlabel("Initial Ash Thickness (cm)")
plt.ylabel("Phosphorus Released (tons)")
plt.title("Phosphorus Release vs Initial Ash Thickness")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("phosphorus_release_vs_thickness_tons.jpg", dpi=300)

# === Fig3===
plt.figure(figsize=(9, 6))
plt.fill_between(T0_cm_array,
                 OC_mean_total - sigma * OC_std_total,
                 OC_mean_total + sigma * OC_std_total,
                 color='lightcoral', alpha=0.5, label='OC Total ±1σ')
plt.plot(T0_cm_array, OC_mean_total, color='darkred', label='OC Total Mean')

plt.xlabel("Initial Ash Thickness (cm)")
plt.ylabel("Organic Carbon Produced (tons)")
plt.title("Organic Carbon Produced vs Initial Ash Thickness (mass ratio)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("organic_carbon_vs_thickness_tons_mass_ratio.jpg", dpi=300)

plt.show()

print("End")
print("- ash_volume_vs_thickness.jpg")
print("- phosphorus_release_vs_thickness_tons.jpg")
print("- organic_carbon_vs_thickness_tons_mass_ratio.jpg")
