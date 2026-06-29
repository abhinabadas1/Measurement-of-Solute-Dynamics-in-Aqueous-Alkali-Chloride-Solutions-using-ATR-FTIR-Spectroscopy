import tkinter as tk
from tkinter import filedialog


def print_world():
    print("Hello World")




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
#air data
air_data=file_read(get_file_path())
#water data
water_1=file_read(get_file_path())
water_2=file_read(get_file_path())
water_3=file_read(get_file_path())

#sample data
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

print(avg_water)
print(avg_sample)

#Now real calculations will begin

