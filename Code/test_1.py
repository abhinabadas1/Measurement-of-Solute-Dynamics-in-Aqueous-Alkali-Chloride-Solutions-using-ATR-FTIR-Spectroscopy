import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#function to convert .txt file data to a 2d array data
def file_read(path):
    # read the file_path file and then I am going to open the file
    try:
        # with is basically a context manager and what open is gonna do is that it is going to open the file which I want to read
        with open(path, 'r', encoding='utf-8') as file:
            data = []
            for line in file:
                #split handles multiple spaces and removes the newline character
                #map iterates the entire the file
                row = list(map(float, line.split()))
                if row:
                    data.append(row)
            return data

    except FileNotFoundError:
        return "Error: File not found"
    except Exception as e:
        # f string is a formatted string, lets you put any variable string into the python code
        return f"an error occurred {e}"

#Now we are going to make it as such that we can select the file path we want using tinker
def get_file_path():
    #create temporary root window
    root=tk.Tk()
    #hide the root window so only the dialog appears
    root.withdraw()
    #Open the openfile Dialog
    file_selected=filedialog.askopenfilename(
        initialdir="/Final Year Project/Code/FIR 26_03_26",
        title="Select your Data file",
        filetypes=(("Text files","*.txt"),("All files", "*.*"))
    )
    #destroy root window after selection
    root.destroy()
    return file_selected


avg_water=[]
avg_sample=[]

print("Select Air Data...")
air_data=file_read(get_file_path())

print("Select Water Data 1, 2, and 3...")
water_1=file_read(get_file_path())
water_2=file_read(get_file_path())
water_3=file_read(get_file_path())

print("Select Sample Data 1, 2, and 3...")
sample_1=file_read(get_file_path())
sample_2=file_read(get_file_path())
sample_3=file_read(get_file_path())

#taking the average of water
for i in range(len(water_1)):
    x=(water_1[i][0]+water_2[i][0]+water_3[i][0])/3
    y=(water_1[i][1]+water_2[i][1]+water_3[i][1])/3
    avg_water.append([x,y])

#taking the average of the sample
for i in range(len(sample_1)):
    x=(sample_1[i][0]+sample_2[i][0]+sample_3[i][0])/3
    y=(sample_1[i][1]+sample_2[i][1]+sample_3[i][1])/3
    avg_sample.append([x,y])

# print(avg_water)
# print(avg_sample)
print("Files loaded and averaged. Starting calculations...")


def get_compound_params(compound_name):

    #Returns the refractive index (n_sample) and initial curve-fitting guesses
    #for specific 1M aqueous alkali chloride solutions.

    #Format of initial_guess: [a1, a2, a3, w1, w2, w3, vd1, vd2, vd3]
    #vd1 = HB Stretching Center
    #vd2 = In-Plane Libration Center
    #vd3 = Out-of-Plane Libration Center


    config = {
        "1M_LiCl": {
            "n_sample": 1.336,
            # Shifts higher due to strong kosmotropic effect (tight water binding)
            "initial_guess": [3000, 9000, 4000, 60, 250, 300, 180, 390, 530]
        },
        "1M_NaCl": {
            "n_sample": 1.337,
            # Standard baseline for moderate hydration
            "initial_guess": [3000, 9000, 4000, 60, 250, 300, 170, 370, 500]
        },
        "1M_KCl": {
            "n_sample": 1.336,
            # Shifts lower due to chaotropic effect (network breaking)
            "initial_guess": [3000, 9000, 4000, 60, 250, 300, 160, 350, 480]
        },
        "1M_CsCl": {
            "n_sample": 1.338,
            # Shifts lowest due to strong chaotropic effect
            "initial_guess": [3000, 9000, 4000, 60, 250, 300, 150, 330, 450]
        }
    }

    if compound_name not in config:
        raise ValueError(f"Compound '{compound_name}' not found. Choose from: {list(config.keys())}")

    return config[compound_name]["n_sample"], config[compound_name]["initial_guess"]

#Now real calculations will begin

# Extracted Constants from Mathematica ---
N_CRYSTAL = 2.38
N_WATER = 1.33

TARGET = "1M_LiCl"

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
    #Normalizes intensity and applies the ATR penetration depth correction
    norm_data = (i_data * 1000) / v
    norm_air = (i_air * 1000) / v
    # dp calculation from your Mathematica file
    dp = (1.0 / v) / ((2 * np.pi) * np.sqrt((N_CRYSTAL**2 / 2.0) - (n_sample**2)))
    return (norm_data - norm_air) / dp

def moving_average(x, w=20):
    #Replicates Mathematica's MovingAverage function
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
    #Single damped harmonic oscillator function.
    numerator = a * w * (v**2)
    term1 = ((v**2) * (w**2)) / (np.pi**2)
    term2 = (vd**2 + (w**2) / (4 * np.pi**2) - v**2)**2
    denominator = 4 * np.pi**2 * (term1 + term2)
    return numerator / denominator

def full_model(v, a1, a2, a3, w1, w2, w3, vd1, vd2, vd3):
    #Combined 3-peak oscillator model
    return (oscillator_peak(v, a1, w1, vd1) +
            oscillator_peak(v, a2, w2, vd2) +
            oscillator_peak(v, a3, w3, vd3))

# Select the specific wavenumber region for fitting (110 to 550 as per the script)
mask = (v_smooth >= 110) & (v_smooth <= 550)
v_fit = v_smooth[mask]
y_fit = delta_alpha_smooth[mask]



bounds = (0, np.inf) # Prevents negative frequencies during optimization

try:
    popt, pcov = curve_fit(full_model, v_fit, y_fit, p0=initial_guess, bounds=bounds, maxfev=10000)
    print("\n--- Curve Fitting Successful ---")
    print(f"HB Stretching Center (vd1): {popt[6]:.2f} cm^-1")
    print(f"In Plane Libration Center (vd2): {popt[7]:.2f} cm^-1")
    print(f"Out of Plane Libration Center (vd3): {popt[8]:.2f} cm^-1")
except Exception as e:
    print(f"\nCurve fitting failed to converge: {e}")
    popt = initial_guess # Fallback to initial guess so the plot still works


# Plotting ---
plt.figure(figsize=(10, 6))

# Plot the smoothed Delta Alpha
plt.plot(v_smooth, delta_alpha_smooth, color='green', label='Delta Alpha (Smoothed)', linewidth=2)

# Plot the fitted model
plt.plot(v_fit, full_model(v_fit, *popt), color='red', linestyle='--', label='Fitted Model')

# Plot the individual deconvoluted peaks
plt.fill_between(v_fit, 0, oscillator_peak(v_fit, popt[0], popt[3], popt[6]), alpha=0.3, label="HB Stretching")
plt.fill_between(v_fit, 0, oscillator_peak(v_fit, popt[1], popt[4], popt[7]), alpha=0.3, label="In Plane Libration")
plt.fill_between(v_fit, 0, oscillator_peak(v_fit, popt[2], popt[5], popt[8]), alpha=0.3, label="Out of Plane Libration")

# Formatting to match standard FTIR plots
plt.xlim(100, 700)
plt.xlabel("Frequency ($cm^{-1}$)", fontsize=12)
plt.ylabel("$\Delta\alpha\ (cm^{-1})$", fontsize=12)
plt.title(f"ATR-FTIR Solute Dynamics (Refractive Index: {N_SAMPLE})", fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()