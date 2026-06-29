import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import glob


def read_ftir_data(file_path):
    # Reads text file containing numeric data and returns it as a list of rows
    data = []
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                # Convert the line into a list of floats, ignoring empty lines
                row = [float(val) for val in line.split()]
                if row:
                    data.append(row)
        return data

    except FileNotFoundError:
        print(f"Warning: File not found - {file_path}")
        return None
    except Exception as e:
        print(f"Warning: Could not read {file_path}. Error: {e}")
        return None


def select_data_folder():
    # Opens a dialog for me to select their data directory
    root = tk.Tk()
    root.withdraw()  # Hide the empty root window

    folder_selected = filedialog.askdirectory(
        initialdir="/Final Year Project/Code/FIR 26_03_26",
        title="Select your FTIR Data Folder"
    )

    root.destroy()
    return folder_selected


# Main script execution

folder_path = select_data_folder()

if folder_path:
    # Grab all the text files in the chosen directory
    search_pattern = os.path.join(folder_path, "*.txt")
    file_list = glob.glob(search_pattern)

    if not file_list:
        print("No .txt files found in the selected folder.")

    for file_path in file_list:
        sample_data = read_ftir_data(file_path)
        file_name = os.path.basename(file_path)

        # Skip to the next file if this one is empty or threw an error
        if not sample_data:
            print(f"Skipping {file_name}: No valid data found.")
            continue

        # Unpack the columns into domain-specific variables
        wavenumbers = [row[0] for row in sample_data]
        transmittance = [row[1] for row in sample_data]

        plt.figure(figsize=(9, 5))
        plt.plot(wavenumbers, transmittance, color='tab:blue', linewidth=1.5, label='Sample Spectrum')

        # Standard practice for FTIR: High to low wavenumbers (left to right)
        plt.gca().invert_xaxis()
        plt.xlim(max(wavenumbers), min(wavenumbers))

        plt.xlabel('Wavenumber (cm⁻¹)', fontsize=11)
        plt.ylabel('Transmittance (%)', fontsize=11)
        plt.title(f'ATR-FTIR Spectrum: {file_name}', fontsize=12)

        plt.grid(True, linestyle=':', alpha=0.7)
        plt.legend()
        plt.tight_layout()

    # rendering all of the plottd spectra
    plt.show()

else:
    print("Folder selection was cancelled. Exiting program.")