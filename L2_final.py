
import tkinter as tk
from tkinter import filedialog
import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from tqdm import tqdm
import concurrent.futures
from PIL import Image, ImageTk
from tkinter import ttk
from tkinter.ttk import Progressbar
from tkinter import Scale
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import multiprocessing
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties
import re
import time

get_ipython().run_line_magic('matplotlib', 'inline')

# Define a global variable to store the paths of extracted text files
extracted_files = []

# Existing code for GUI setup, function definitions, etc.

def choose_folder_l2():
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    folder_path = filedialog.askdirectory()  # Open a dialog to select a folder
    return folder_path
 
def plot_selected_files_l2():
    # Create a loading message window
    loading_message = tk.Toplevel(root)
    loading_message.title("Loading...")
    loading_message.geometry("200x100")
    loading_label = tk.Label(loading_message, text="Loading, please wait...")
    loading_label.pack()
    
    # Force the window to update and display
    loading_message.update()
    
    # Get the selected files from the listbox
    selected_file_indices = listbox.curselection()
    
    if selected_file_indices:
        # Clear the plot before plotting the new files
        plt.clf()
        
        print('file selected')
        
        # Initialize arrays to store longitude, latitude,time and selected variables
        merged_lons = []
        merged_lats = []
        merged_datetime = []
        merged_variable = []
        
        def plot_file(selected_file_info):
            
            nonlocal merged_lons, merged_lats, merged_datetime, merged_variable
            
            selected_file_path = selected_file_info[0]
            print("Plotting selected file:", selected_file_path)
            
            print(selected_variable)
               
            # Read the NetCDF file using xarray
            ds_expert = xr.open_dataset(selected_file_path)
            ds_expert = ds_expert.assign_coords(longitude=(((ds_expert.longitude + 180) % 360) - 180))
            
            
            # Extract longitude, latitude,time and variable values
            try:
                print(selected_variable)
                lon = ds_expert['longitude'].values
                
                lat = ds_expert['latitude'].values
                
                time = ds_expert['time'].values
                
                variable_data = ds_expert[selected_variable].values
                
            except Exception as e:
                
                print("Exiting -- Error occured")
                print(e)
 
            # Append longitude, latitude,time and selected variables to the merged arrays
            merged_lons.extend(lon)
            merged_lats.extend(lat)
            merged_datetime.extend(time)
            merged_variable.extend(variable_data)
            
        # Get selected file info for all selected files
        selected_files_info = (listbox.get(index) for index in selected_file_indices)
        
        selected_variable = variable_dropdown.get()  # Get the selected variable from the dropdown
        
        print(selected_variable)
        
        # Plot selected files in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor() as executor:
            
            executor.map(plot_file, selected_files_info,timeout=5)
     
         
        # Convert merged lists to numpy arrays
        merged_lons = np.array(merged_lons)
        merged_lats = np.array(merged_lats)
        merged_datetime = np.array(merged_datetime)
        merged_variable = np.array(merged_variable)
                 
        # Get latitude and longitude and time ranges from entry fields
        min_lat = float(min_lat_entry.get())
        max_lat = float(max_lat_entry.get())
        min_lon = float(min_lon_entry.get())
        max_lon = float(max_lon_entry.get())

        
        start_date = str(start_date_entry.get())
        end_date = str(end_date_entry.get())
        
        start_time = datetime.strptime(start_date, '%Y-%m-%d')
        end_time = datetime.strptime(end_date, '%Y-%m-%d')
        
        min_variable, max_variable = float(min_variable_entry.get()), float(max_variable_entry.get())
        
        if selected_variable:
            print('plot..')
            merged_fig = plot_data_l2(merged_lons, merged_lats, merged_variable, min_lat, max_lat,min_lon,max_lon, min_variable,max_variable)
            update_plot_canvas_l2(merged_fig)   # Update the plot on the canvas
            
            # Close the loading message window after plotting is completed
            loading_message.destroy()
        
    else:
        print("No file selected.....")
      
selected_files = []
latitudes = []
longitudes = []
date_files = []

def filter_files_by_coordinates_l2(folder_path, min_lat, max_lat,min_lon,max_lon,start_time,end_time,min_variable,max_variable,selected_variable, pbar=None):
    #for calculating time
    start_time_filter_files = time.time()
        
    #code for the function 
    valid_files_list = []
    time_list = []
    
    # List to store filtered file paths
    filtered_files = []
    
    # Total number of files
    file_list = os.listdir(folder_path)
    all_files = len(os.listdir(folder_path))
        
    # Iterate through files in the folder
    for file_name in file_list:
        # Define the regular expression pattern to match the date in the file name
        
        pattern = re.compile(r'(\d{8})')
        
        # Find all matches of the pattern in the file name
        matches = re.findall(pattern, file_name)
        #fde=datetime(2023,10,10)
        if matches:
            date_str = matches[0]
            
            # Format the date as YYYY-MM-DD
            date = '-'.join([date_str[:4], date_str[4:6], date_str[6:]])
            
            file_date = datetime.strptime(date, '%Y-%m-%d')
            
            #comparing the time that is between the user given start time and end time
            if start_time <=  file_date  <= end_time:
                selected_files.append(os.path.join(folder_path, file_name))
                print(file_date)

    file_length  = len(selected_files)
               
    #iterate through selected files
    for idx,file_name in enumerate(selected_files,1):
        
        # Update progress bar
        if pbar is not None:
            pbar['value'] = idx / file_length * 100
            pbar.update()
            
        # Check if the file is a NetCDF file
        if file_name.endswith('.nc'):
            
            
            try:
                
                # Open the NetCDF file
                with xr.open_dataset(file_name) as dataset:
                    
                    # Extract latitude and longitude,time and variable coordinates from the dataset
                    
                    latitudes = dataset['latitude'].values
                    longitudes = dataset['longitude'].values
                    nc_time = dataset['time'].values
                    variable_data = dataset[selected_variable].values
                    print("Variable data : ",variable_data)
                       
                    # Find the indices corresponding to latitude and longitude and time ranges
                    lat_indices = np.where((latitudes >= min_lat) & (latitudes <= max_lat))[0]
                    lon_indices = np.where((longitudes >= min_lon) & (longitudes <= max_lon))[0]
                    variable_indices = np.where((variable_data >= min_variable) & (variable_data <= max_variable))[0]
                    

                    #Check for common indices (intersection) for latitude and longitude 
                    common_indices = np.intersect1d(np.intersect1d(lat_indices, lon_indices),variable_indices)
                    
                    
                    if len(common_indices) != 0:

                        filtered_latitudes = latitudes[common_indices]
                        filtered_longitudes = longitudes[common_indices]
                        filtered_variable = variable_data[common_indices]
                        
                        # Check if there are non-NaN variable values within the specified latitude and longitude and time ranges
                        if not np.isnan(filtered_variable).all():
                            filtered_files.append(( file_name, filtered_latitudes, filtered_longitudes,filtered_variable))
                        
            except(OSError, RuntimeError) as e:
                print(f"Error processing file : {e}")
                continue
            
    end_time_filter_files = time.time()
    elapsed_time = end_time_filter_files - start_time_filter_files
    print(f"The filter_files function took {elapsed_time:.2f} seconds to run.")
            
    return filtered_files


def button_3_clicked_l2():
    
    # calculating time
    start_time_button_3 = time.time()
    
    # code for the function
    
    # Get latitude and longitude values and time values from the entry fields
    min_lat = float(min_lat_entry.get())
    max_lat = float(max_lat_entry.get())
    min_lon = float(min_lon_entry.get())
    max_lon = float(max_lon_entry.get())
    
    start_date = str(start_date_entry.get())
    end_date = str(end_date_entry.get())
    
    start_time = datetime.strptime(start_date, '%Y-%m-%d')
    end_time = datetime.strptime(end_date, '%Y-%m-%d')
    
    min_variable, max_variable = float(min_variable_entry.get()), float(max_variable_entry.get())
        
    # Filter files based on the provided coordinates
    folder_path = choose_folder_l2()  # Choose the folder path through file management system
    
    selected_variable = variable_dropdown.get()  # Get the selected variable from the dropdown
    
    # Remove previous progress bar if it exists
    for widget in root.winfo_children():
        if isinstance(widget, Progressbar):
            widget.destroy()
            
    # create a new progress bar
    pbar = Progressbar(root, orient="horizontal", length=200, mode="determinate")
    pbar.pack(side="bottom")
    
    filtered_files = filter_files_by_coordinates_l2(folder_path, min_lat,max_lat,min_lon,max_lon, start_time,end_time,min_variable, max_variable, selected_variable, pbar)
       
    valid_files = len(filtered_files)
    
    print('Filtered Files: ')
    
    for file_name in filtered_files:
        print(file_name)
    print('number of valid files in the given lat,lon and time range : ', valid_files)
    print("Total files:", len(os.listdir(folder_path)))
    
    # Update the listbox with the filtered files
    listbox.delete(0, tk.END)  # Clear the listbox
    for file_name in filtered_files :
        listbox.insert(tk.END, file_name)
        
    # Update the label with the count of valid files
    valid_files_label.config(text=f"Valid Files: {valid_files}/{len(os.listdir(folder_path))}")
    
    end_time_button_3 = time.time()
    elapsed_time = end_time_button_3 - start_time_button_3
    print(f"The button_3  function took {elapsed_time:.2f} seconds to run.")
    
def plot_data_l2(lons, lats, variable_data,min_lat, max_lat, min_lon, max_lon, min_variable,max_variable):
    #calculating time
    start_time_plot = time.time()
    
    #code for the function
    
    # Define the path to the coastline shapefile
    coastline_shapefile = "/home/guest/Debangshu_SWOT/Debangshu_SWOT_GUI/ne_10m_coastline.shp"
    
    # Check if the shapefile exists
    if os.path.exists(coastline_shapefile):
        # Read the coastline shapefile
        coastlines = gpd.read_file(coastline_shapefile)
        
        # Create a figure and axis using Cartopy
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw=dict(projection=ccrs.PlateCarree()))
         
        # Add blue marble background
        ax.stock_img()
        
        # Plot variable values within the specified latitude and longitude and time ranges
        pcm = ax.scatter(lons, lats, c=variable_data, cmap='Spectral_r', vmin=min_variable, s = 10 , vmax=max_variable, transform=ccrs.PlateCarree())
        ax.set_extent([min_lon,max_lon,min_lat,max_lat], crs=ccrs.PlateCarree())
        
        # add coastlines
        coastlines.plot(ax=ax, color='black', linewidth=1.5)
        
        # Add additional features like land, borders, etc.
        #ax.add_feature(cfeature.LAND, facecolor='lightgray')
        #ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='black')
        #ax.add_feature(cfeature.LAKES, edgecolor='blue', facecolor='blue', alpha=0.5)
        #ax.add_feature(cfeature.RIVERS, edgecolor='blue')

           
        # add colorbar
        cbar = plt.colorbar(pcm, ax=ax, orientation='vertical', shrink=0.7, label='Selected Variable (m)')
        
        # add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        
       
        # Add title and labels
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    
        
        # return the fig
        return fig
    else:
        print("Coastline shapefile not found.")
        
    end_time_plot = time.time()
    elapsed_time = end_time_plot - start_time_plot
    print(f"The plot function took {elapsed_time:.2f} seconds to run.")
    
        
def update_plot_canvas_l2(fig):
    
    global canvas_new
    # Clear the previous plot on canvas
    for widget in center_frame.winfo_children():
        widget.destroy()
        
    # Create a new matplotlib canvas for embedding the updated plot
    canvas_new = FigureCanvasTkAgg(fig, master=center_frame)
    canvas_new.draw()
    canvas_new.get_tk_widget().pack(expand=True, fill=tk.BOTH)
    
    # add matplotlib navigation toolbar
    toolbar = NavigationToolbar2Tk(canvas_new, center_frame)
    toolbar.update()
    canvas_new._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    
def extract_data_to_text_l2(nc_file_path, output_folder, min_lat,max_lat,min_lon,max_lon, start_time,end_time,min_variable,max_variable):
    try:
        # Open the NetCDF file using xarray
        ds = xr.open_dataset(nc_file_path)
        selected_variable = variable_dropdown.get()  # Get the selected variable from the dropdown
        
        # Extract latitude, longitude, and SSHA data
        latitudes = ds['latitude'].values
        longitudes = ds['longitude'].values
        plot_variable_values = ds[selected_variable].values
        
        # Extract time data
        time_data = ds['time'].values.astype('datetime64[s]').astype(datetime)  # Convert to datetime object
        
        # Flatten arrays and get corresponding indices of non-NaN values in SSHA
        valid_indices = np.where(~np.isnan(plot_variable_values))

        # Filter data based on user-defined latitude and longitude ranges
        latitudes = latitudes[valid_indices].flatten()  # Ensure 1-dimensional
        longitudes = longitudes[valid_indices].flatten()  # Ensure 1-dimensional
        #time_data = time_data[valid_indices].flatten()
        plot_variable_values = plot_variable_values[valid_indices].flatten()  # Ensure 1-dimensional

        latitudes = np.round(latitudes, 2)  # Round latitude values to two decimal places
        longitudes = np.round(longitudes, 2)  # Round longitude values to two decimal places
        plot_variable_values = np.round(plot_variable_values, 2)  # Round SSHA values to two decimal places
        
        
        # Convert time data to the desired format
        formatted_time_data = []
        for time in time_data:
            formatted_time = time.strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]  # Convert datetime to string with milliseconds
            formatted_time = formatted_time.replace('-', '   ').replace(':', '   ').replace('.', '   ')  # Replace '-' ':' with spaces
            formatted_time = '  '.join(formatted_time.split())  # Remove extra spaces
            formatted_time_data.append(formatted_time)
            
        text_date_format= []
            
        for date_string in formatted_time_data:
            # Parse the date string
            date_obj = datetime.strptime(date_string, "%Y %m %d %H %M %S")
            
            # Format the date as YYYY-MM-DD
            formatted_date = date_obj.strftime("%Y-%m-%d")
            
            text_date_format.append(formatted_date)
  
        print(f"Time data formatted: {formatted_time_data}")
    
        # Construct the output file path
        output_file_path = os.path.join(output_folder, os.path.basename(nc_file_path).replace('.nc', '_extracted_data.txt'))
        
        # Write latitude, longitude, SSHA, and time data to the text file
        with open(output_file_path, 'w') as file:
            
            for time, lat, lon, variable in zip(text_date_format, latitudes, longitudes, plot_variable_values):
                file.write(f"{time}\t{lat}\t{lon}\t{variable}\n")

        print(f"Extracted data from {nc_file_path} saved to: {output_file_path}")

        # Return the file path
        return output_file_path
    
    except Exception as e:
        print(f"Error extracting data from {nc_file_path}: {e}")
        return ''  # Return an empty string if an error occurs


def save_selected_files_as_text_l2():
    # Get the indices of the selected files in the listbox
    selected_indices = listbox.curselection()
    
    if selected_indices:
        # Create a folder to save the text files if it doesn't exist
        output_folder = filedialog.askdirectory(title="Select output folder")
        print(f"Output folder: {output_folder}")  # Print output folder path
        
        if output_folder:
            # Prompt the user for the custom file name
            custom_file_name = tk.simpledialog.askstring("Custom File Name", "Enter a custom file name:")
            if custom_file_name:
                combined_file_name = f"{custom_file_name}.txt"
            else:
                combined_file_name = "combined_files.txt"
                
            # Initialize a list to store file paths of all extracted text files
            all_files_paths = []
            
            min_lat = float(min_lat_entry.get())
            max_lat = float(max_lat_entry.get())
            min_lon = float(min_lon_entry.get())
            max_lon = float(max_lon_entry.get())
            
            start_date = str(start_date_entry.get())
            end_date = str(end_date_entry.get())
            
            start_time = datetime.strptime(start_date, '%Y-%m-%d')
            end_time = datetime.strptime(end_date, '%Y-%m-%d')
            
            min_variable, max_variable = float(min_variable_entry.get()), float(max_variable_entry.get())
                   
            # Initialize an empty string to hold the combined content of all text files
            combined_content = ""
            
            for index in selected_indices:
                file_info = listbox.get(index)
                nc_file_path = file_info[0]
                time_array = []  # Placeholder for time array 
                file_path = extract_data_to_text_l2(nc_file_path, output_folder, min_lat,max_lat,min_lon,max_lon, start_time,end_time, min_variable, max_variable)
                if file_path:
                    all_files_paths.append(file_path)
                    
            # Combine content of all extracted text files
            for file_path in all_files_paths:
                with open(file_path, 'r') as file:
                    file_content = file.read()  # Read content of each file
                    combined_content += file_content + '\n'  # Add a newline for separation between files
                    print(f"Content of {file_path}: {file_content}")  # Print content of each file
           
            # Write the combined content to a single file
            combined_file_path = os.path.join(output_folder, combined_file_name)
            with open(combined_file_path, "w") as combined_file:
                combined_file.write(combined_content)
           
            print(f"Text files saved and combined successfully as: {combined_file_path}")
    else:
        print("No file selected.")
        

# Create the main window
root = tk.Tk()
root.title("SWOT_GUI_level 2")
root.minsize(1660, 1000)

def set_default_l2(entry, value):
    entry.delete(0, tk.END)  # Clear any existing text
    entry.insert(0, value)    # Insert default date

# Load and resize the first image (adjust path and size as needed)
image_path1 = '/home/guest/Debangshu_SWOT/Debangshu_SWOT_GUI/download.jpeg'
image1 = Image.open(image_path1).resize((60, 60))  # Adjust the resize dimensions
photo1 = ImageTk.PhotoImage(image1)

# Load and resize the second image (adjust path and size as needed)
image_path2 = "/home/guest/Debangshu_SWOT/Debangshu_SWOT_GUI/Indian_Space_Research_Organisation_Logo.svg.png"
image2 = Image.open(image_path2).resize((60, 60))  # Adjust the resize dimensions
photo2 = ImageTk.PhotoImage(image2)

# Create top frame with sky blue background
top_frame = tk.Frame(root, width=200, height=40, pady=10, bg="sky blue")
top_frame.pack(side="top", fill="both")

# Add the first image to the top-left corner of the top frame
image_label1 = tk.Label(top_frame, image=photo1)
image_label1.image = photo1  # Store a reference to the image
image_label1.pack(side="left", anchor="nw")

# Add the second image next to the first image in the top frame
image_label2 = tk.Label(top_frame, image=photo2)
image_label2.image = photo2  # Store a reference to the image
image_label2.pack(side="left", anchor="nw")

# Load and display the GIF image in the left frame
gif_path = "/home/guest/Debangshu_SWOT/Debangshu_SWOT_GUI/karin_animation.gif"  # Update the path to your GIF image
gif_image = Image.open(gif_path)


# Resize the GIF image to a smaller dimension
new_width = 50  # Adjust the width as needed
new_height = 50  # Adjust the height as needed
gif_resized = gif_image.resize((new_width, new_height), Image.BICUBIC)  # Use Image.BICUBIC as an alternative


# Convert the resized image to a Tkinter PhotoImage
gif_photo = ImageTk.PhotoImage(gif_resized)

# Calculate the width of the left frame
left_frame_width = new_width + 60  # Add extra padding as needed

# Create the left frame with the adjusted width
left_frame = tk.Frame(root, width=left_frame_width, height=new_height + 60, padx=30, bg="sky blue")
left_frame.pack(side="right", fill="both")

# Create a label to display the resized GIF image
gif_label = tk.Label(left_frame, image=gif_photo, bg="sky blue")
gif_label.pack(side="top", pady=20)


def update_image_l2(ind):
    frame = gif_image.seek(ind)
    gif_frame = ImageTk.PhotoImage(gif_image)
    gif_label.configure(image=gif_frame)
    gif_label.image = gif_frame
    ind += 1
    root.after(100, update_image_l2, ind % gif_image.n_frames)
    

# Start the GIF animation
update_image_l2(0)

# Create fields for latitude and longitude and time with bold font
font_bold = ("Arial", 10, "bold")

# First line: Latitude labels and entries
latitude_frame = tk.Frame(left_frame, bg="sky blue")
latitude_frame.pack(side="top", fill="x", padx=5, pady=5)

min_lat_label = tk.Label(latitude_frame, text="Min Latitude:", font=font_bold)
min_lat_default = 12.0 # set default minimum latitude value
min_lat_entry = tk.Entry(latitude_frame, font=font_bold)
set_default_l2(min_lat_entry, min_lat_default)


max_lat_label = tk.Label(latitude_frame, text="Max Latitude:", font=font_bold)
max_lat_default = 30.0 # set default maximum latitude value
max_lat_entry = tk.Entry(latitude_frame, font=font_bold)
set_default_l2(max_lat_entry, max_lat_default)

min_lat_label.pack(side="left", padx=45, pady=5)
min_lat_entry.pack(side="left", padx=5, pady=5)
max_lat_label.pack(side="left", padx=10, pady=5)
max_lat_entry.pack(side="left", padx=5, pady=5)


# Second line: Longitude labels and entries
longitude_frame = tk.Frame(left_frame, bg="sky blue")
longitude_frame.pack(side="top", fill="x", padx=5, pady=5)

min_lon_label = tk.Label(longitude_frame, text="Min Longitude:", font=font_bold)
min_lon_default = 50.0 # set default minimum longitude value
min_lon_entry = tk.Entry(longitude_frame, font=font_bold)
set_default_l2(min_lon_entry, min_lon_default)

max_lon_label = tk.Label(longitude_frame, text="Max Longitude:", font=font_bold)
max_lon_default = 75.0 # set default maximum longitude value
max_lon_entry = tk.Entry(longitude_frame, font=font_bold)
set_default_l2(max_lon_entry, max_lon_default)

min_lon_label.pack(side="left", padx=40, pady=5)
min_lon_entry.pack(side="left", padx=5, pady=5)
max_lon_label.pack(side="left", padx=10, pady=5)
max_lon_entry.pack(side="left", padx=5, pady=5)

# third line : date entries   
time_frame = tk.Frame(left_frame, bg = 'sky blue')
time_frame.pack(side='top',fill='x',padx=5,pady=5)

start_date_label = tk.Label(time_frame, text = 'Start date: ',font = font_bold)
default_start_date = '2023-07-26' # set default start date
start_date_entry = tk.Entry(time_frame, font = font_bold)
set_default_l2(start_date_entry, default_start_date)


end_date_label = tk.Label(time_frame, text = 'End Date: ', font = font_bold)
default_end_date = '2023-07-27' # set default end date
end_date_entry = tk.Entry(time_frame,font = font_bold)
set_default_l2(end_date_entry, default_end_date)

start_date_label.pack(side ='left',padx=50,pady=5)
start_date_entry.pack(side = 'left',padx=5,pady=5)
end_date_label.pack(side='left',padx=15,pady=5)
end_date_entry.pack(side='left',padx=5,pady=5)

# fourth line : selected variable min max
variable_frame = tk.Frame(left_frame,bg = 'sky blue')
variable_frame.pack(side = 'top',fill='x',padx = 5,pady = 5)

min_variable_label = tk.Label(variable_frame,text = 'Variable min',font = font_bold)
min_variable_default = -0.5 # set default minimum variable value
min_variable_entry = tk.Entry(variable_frame,font = font_bold)
set_default_l2(min_variable_entry, min_variable_default)

max_variable_label = tk.Label(variable_frame,text = 'Variable max',font = font_bold)
max_variable_default = 0.5 # set default maximum variable value
max_variable_entry = tk.Entry(variable_frame,font = font_bold)
set_default_l2(max_variable_entry, max_variable_default)

min_variable_label.pack(side='left',padx = 40,pady = 5)
min_variable_entry.pack(side='left',padx = 5,pady = 5)
max_variable_label.pack(side='left',padx=10,pady = 5)
max_variable_entry.pack(side='left',padx = 5,pady = 5)


# Centering latitude and longitude ssha and time frames
latitude_frame.pack_configure(anchor="center")
longitude_frame.pack_configure(anchor="center")
time_frame.pack_configure(anchor='center')
variable_frame.pack_configure(anchor='center')

# Create buttons in left frame with bold font and sky blue color
button_font = ("Arial", 10, "bold")
button_bg_color = "sky blue"

# Create a frame to contain the buttons
button_frame_l2 = tk.Frame(left_frame, bg="sky blue")
button_frame_l2.pack(side="top", padx=15, pady=25)

button3_l2 = tk.Button(button_frame_l2, text="Select folder for the data", command=button_3_clicked_l2, font=button_font, bg=button_bg_color)
button3_l2.grid(row=0, column=2, padx=(10, 30), pady=5)  # Increased padding on the right

# Create a frame to contain the listbox
list_frame = tk.Frame(root, bg="white", width=200, height=200)
list_frame.pack(side="bottom", fill="both")

# Create a scrollbar for the listbox
scrollbar = tk.Scrollbar(list_frame, orient="vertical")

# Create a listbox widget to display the filtered files
listbox = tk.Listbox(list_frame, yscrollcommand=scrollbar.set, width=50, height=10, selectmode=tk.MULTIPLE)
scrollbar.config(command=listbox.yview)
scrollbar.pack(side="right", fill="y")
listbox.pack(side="left", fill="both", expand=True)

# Create a label to display the count of valid files
valid_files_label = tk.Label(root, text="Valid Files: 0/0")
valid_files_label.pack(side="bottom")


# variable drop down menu geometry
variable_dropdown_label = tk.Label(button_frame_l2,text = 'Select Variable',font=button_font, bg=button_bg_color)
variable_dropdown_label.grid(row=0, column=0, padx=(20, 10), pady=5)
variable_dropdown = tk.StringVar()
variable_dropdown = ttk.Combobox(button_frame_l2, values=["ssha_karin", "ssha_karin_2","ssh_karin_2",'distance_to_coast','heading_to_coast','geoid'])
variable_dropdown.set("ssha_karin") # setting of the default value
variable_dropdown.grid(row=0, column=1)

# Plot button
plot_variable = tk.Button(button_frame_l2, text="Plot the selected variable", command=plot_selected_files_l2, font=button_font, bg=button_bg_color)
plot_variable.grid(row=0, column=3, padx=(5, 200), pady=5)  # Increased padding on the left

# Create a frame to contain the "Save Selected as Text" button
save_button_l2_frame = tk.Frame(left_frame, bg="sky blue")
save_button_l2_frame.pack(side="top", padx=15, pady=25)

save_text_button_l2 = tk.Button(save_button_l2_frame, text="Save Selected as Text", command=save_selected_files_as_text_l2, font=button_font, bg=button_bg_color)
save_text_button_l2.grid(row=0, column=0, padx=(30, 10), pady=0)

# Create center frame for displaying plots
center_frame = tk.Frame(root)
center_frame.pack(expand=True, fill="both")

# Start the GUI event loop
root.mainloop()
