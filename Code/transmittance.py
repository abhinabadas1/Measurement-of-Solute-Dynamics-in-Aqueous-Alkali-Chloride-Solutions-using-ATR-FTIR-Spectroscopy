import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt


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

sample=file_read(get_file_path())

x = [row[0] for row in sample]
y = [row[1] for row in sample]

plt.figure(figsize=(9, 5))

# Plot the spectrum with standard Matplotlib defaults
plt.plot(x, y, color='tab:blue', linewidth=1.5, label='Sample Spectrum')

# Reverse the X-axis (Standard practice for all FTIR data)
plt.gca().invert_xaxis()

# Set axis limits strictly to the data bounds
plt.xlim(max(x), min(x))

# Labels and Title
plt.xlabel('Wavenumber (cm⁻¹)', fontsize=11)
plt.ylabel('Transmittance (%)', fontsize=11)
plt.title('ATR-FTIR Spectrum', fontsize=12)

# Standard grid and legend
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()

plt.tight_layout()
plt.show()

