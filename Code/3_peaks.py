import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Function to convert .txt file data to a 2d array data
def file_read(path):
    try:
        with open(path, 'r', encoding='utf-8') as file:
            data = []
            for line in file:
                row = list(map(float, line.split()))
                if row:
                    data.append(row)
            return data
    except FileNotFoundError:
        return "Error: File not found"
    except Exception as e:
        return f"an error occurred {e}"

def get_file_path():
    root = tk.Tk()
    root.withdraw()
    file_selected = filedialog.askopenfilename(
        initialdir="/Final Year Project/Code/FIR 26_03_26",
        title="Select your Data file",
        filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
    )
    root.destroy()
    return file_selected


avg_water = []
avg_sample = []

print("Select Air Data...")
air_data = file_read(get_file_path())

print("Select Water Data 1, 2, and 3...")
water_1 = file_read(get_file_path())
water_2 = file_read(get_file_path())
water_3 = file_read(get_file_path())

print("Select Sample Data 1, 2, and 3...")
sample_1 = file_read(get_file_path())
sample_2 = file_read(get_file_path())
sample_3 = file_read(get_file_path())

# Taking the average of water
for i in range(len(water_1)):
    x = (water_1[i][0] + water_2[i][0] + water_3[i][0]) / 3
    y = (water_1[i][1] + water_2[i][1] + water_3[i][1]) / 3
    avg_water.append([x, y])

# Taking the average of the sample
for i in range(len(sample_1)):
    x = (sample_1[i][0] + sample_2[i][0] + sample_3[i][0]) / 3
    y = (sample_1[i][1] + sample_2[i][1] + sample_3[i][1]) / 3
    avg_sample.append([x, y])

print("Files loaded and averaged. Starting calculations...")


def get_compound_params(compound_name):
    """
    Format of initial_guess: [a1, a2, a3, w1, w2, w3, vd1, vd2, vd3]
    a1 is negative to account for the bleach effect in HB stretching.
    """
    config = {
        "1M_LiCl": {
            "n_sample": 1.336,
            "initial_guess": [-3000, 9000, 4000, 60, 250, 300, 180, 390, 530]
        },
        "1M_NaCl": {
            "n_sample": 1.337,
            "initial_guess": [-3000, 9000, 4000, 60, 250, 300, 170, 370, 500]
        },
        "1M_KCl": {
            "n_sample": 1.336,
            "initial_guess": [-3000, 9000, 4000, 60, 250, 300, 160, 350, 480]
        },
        "1M_CsCl": {
            "n_sample": 1.338,
            "initial_guess": [-3000, 9000, 4000, 60, 250, 300, 150, 330, 450]
        }
    }

    if compound_name not in config:
        raise ValueError(f"Compound '{compound_name}' not found. Choose from: {list(config.keys())}")

    return config[compound_name]["n_sample"], config[compound_name]["initial_guess"]


# Extracted Constants from Mathematica ---
N_CRYSTAL = 2.38
N_WATER = 1.33
TARGET = "1M_CsCl" # Update target as necessary

# Fetch the specific parameters for this run
N_SAMPLE, initial_guess = get_compound_params(TARGET)

print(f"Processing {TARGET} with Refractive Index: {N_SAMPLE}")

# Converting the arrays to NumPy arrays for math ---
air_arr = np.array(air_data)
water_arr = np.array(avg_water)
sample_arr = np.array(avg_sample)

# Extract Wavenumber (v) and Intensities
v = air_arr[:, 0]
i_air = air_arr[:, 1]
i_water = water_arr[:, 1]
i_sample = sample_arr[:, 1]


# Optical Constants & Subtraction Functions ---
def calculate_alpha(v, i_data, i_air, n_sample):
    norm_data = (i_data * 1000) / v
    norm_air = (i_air * 1000) / v
    dp = (1.0 / v) / ((2 * np.pi) * np.sqrt((N_CRYSTAL**2 / 2.0) - (n_sample**2)))
    return (norm_data - norm_air) / dp

def moving_average(x, w=20):
    return np.convolve(x, np.ones(w), 'valid') / w

# Execute the physics math
alpha_water = calculate_alpha(v, i_water, i_air, N_WATER)
alpha_sample = calculate_alpha(v, i_sample, i_air, N_SAMPLE)
delta_alpha = alpha_sample - alpha_water

# Apply 20-point moving average
v_smooth = v[19:]
delta_alpha_smooth = moving_average(delta_alpha, w=20)


# Curve Fitting Models (Damped Harmonic Oscillator) ---
def oscillator_peak(v, a, w, vd):
    numerator = a * w * (v**2)
    term1 = ((v**2) * (w**2)) / (np.pi**2)
    term2 = (vd**2 + (w**2) / (4 * np.pi**2) - v**2)**2
    denominator = 4 * np.pi**2 * (term1 + term2)
    return numerator / denominator

def full_model(v, a1, a2, a3, w1, w2, w3, vd1, vd2, vd3):
    return (oscillator_peak(v, a1, w1, vd1) +
            oscillator_peak(v, a2, w2, vd2) +
            oscillator_peak(v, a3, w3, vd3))

# Select the extended wavenumber region for fitting (100 to 600 cm-1)
mask = (v_smooth >= 100) & (v_smooth <= 600)
v_fit = v_smooth[mask]
y_fit = delta_alpha_smooth[mask]

# --- 3-PEAK BOUNDS ---
# a1, a2, a3 (Amplitudes): -infinity to infinity
# w1, w2, w3 (Widths): restricted between 10 and 600
# vd1 (HB Stretch): 120 to 220
# vd2 (In-Plane Libration): 280 to 450
# vd3 (Out-of-Plane Libration): 450 to 600
lower_bounds = [-np.inf, -np.inf, -np.inf, 10, 50, 50, 120, 280, 450]
upper_bounds = [np.inf, np.inf, np.inf, 400, 600, 600, 220, 450, 600]
bounds = (lower_bounds, upper_bounds)

try:
    popt, pcov = curve_fit(full_model, v_fit, y_fit, p0=initial_guess, bounds=bounds, maxfev=10000)
    print("\n--- Curve Fitting Successful ---")
    print(f"HB Stretching Center (Peak 1) (vd1): {popt[6]:.2f} cm^-1")
    print(f"Rattling Mode Center(Peak 2) (vd2): {popt[7]:.2f} cm^-1")
    print(f"Libration Mode Center (Peak 3) (vd3): {popt[8]:.2f} cm^-1")
except Exception as e:
    print(f"\nCurve fitting failed to converge: {e}")
    popt = initial_guess


# Plotting ---
plt.figure(figsize=(10, 6))

# Plot the smoothed Delta Alpha (dotted line)
plt.plot(v_smooth, delta_alpha_smooth, color='green', linestyle=':', label='Delta Alpha (Smoothed)', linewidth=3)

# Plot the fitted model (solid red line)
plt.plot(v_fit, full_model(v_fit, *popt), color='red', linestyle='-', label='Fitted Model')

# Plot the individual deconvoluted peaks for the 3-peak model with custom colors
plt.fill_between(v_fit, 0, oscillator_peak(v_fit, popt[0], popt[3], popt[6]), color='orchid', alpha=0.5, label="Peak 1 (HB Stretching)")
plt.fill_between(v_fit, 0, oscillator_peak(v_fit, popt[1], popt[4], popt[7]), color='lightgreen', alpha=0.5, label="Peak 2 (Rattling Mode)")
plt.fill_between(v_fit, 0, oscillator_peak(v_fit, popt[2], popt[5], popt[8]), color='olive', alpha=0.5, label="Peak 3 (Libration Mode)")

# Formatting adjusted for 600 cm-1 upper limit
plt.xlim(100, 600)
plt.xlabel("Frequency (cm$^{-1}$)", fontsize=12)
plt.ylabel("$\Delta\\alpha$ (cm$^{-1}$)", fontsize=12)
plt.title(f"ATR-FTIR Solute Dynamics of {TARGET} (Refractive Index: {N_SAMPLE})", fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)

plt.show()