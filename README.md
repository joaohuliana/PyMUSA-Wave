# Attention: Working on Python version 3.12.3

# pyMUSA-WAVE: Ultrasound Data Processing Platform

Description:
A software platform for processing ultrasound data, specifically designed for Acoustic Radiation Force Imaging (ARFI)
and Magneto-Motive Ultrasound (MMUS). It provides tools for RF data manipulation, particle displacement/velocity estimation, 
and Shear Wave (SW) phase and group velocity analysis.

Reference:
Uliana, J. H., Sampaio, D. R. T., Sanches, A. F., Pavan, T. Z., & Carneiro, A. A. D. O. (2024).
MUSA-WAVE: A Software Platform for Magnetomotive Ultrasound and Shear Wave Analysis.
2024 IEEE UFFC Latin America Ultrasonics Symposium (LAUS) Proceedings, 7–10. https://doi.org/10.1109/LAUS60931.2024.10553040


# How to Install

1.    Install Python 3.12.3 from https://www.python.org/ftp/python/3.12.3/python-3.12.3-amd64.exe.

2.    Install VSCode (or your preferred IDE) from https://code.visualstudio.com/download.

3.    Download the project (Download ZIP) and extract it.

4.    In VSCode, go to File > Open Folder and select the pyMUSA-Wave-main folder.

5.    Click Terminal > New Terminal and create a virtual environment: python -m venv .venv

6.    Activate the virtual environment: \.venv\Scripts\activate

7.    Install the required modules: pip install -r requirements.txt

8.    Run the pyvma_ppg_white.py file.


# How to Use (Example File)

1.    Download an example file from https://drive.google.com/drive/folders/14WjUnAMyzCzxjlcRgJcm4UGMGeQRocKI.

2.    Go to File > Load .mat file and select CIRS_049_Type_I_1.mat.

        This example contains Acoustic Radiation Force Imaging (ARFI) data from a CIRS 049 phantom using an ATL L4-7 transducer.

        The default settings in the "US Parameters" section match the example.

 3.   Processing RF Data (Left Column): Controls for processing ultrasound RF data.

 4.   Generating Displacement Maps (Middle Column): Generates and processes displacement maps.

 5.   Click the "Loupas" button and wait for the process to finish.

 6.   The "particle velocity map" can be navigated frame by frame using the horizontal slider.

 7.   The vertical sliders control the colorbar's minimum and maximum values.

 8.   Click "Low-Pass" to remove high-frequency noise, then click "Lee's" to despeckle the image.

 9.   Shear Wave Analysis (Right Column): To calculate the Shear Wave group velocity, click on the "2D map" in the right column.

 10.   As with the displacement maps, the vertical sliders control the colorbar.

 11.   The "cursor" button displays the value at each point on the SW group velocity map.

 12.   The combobox offers "Vx" (horizontal), "Vy" (vertical), and "Velocity" (magnitude) options for SW velocity values.

 13.   The "ROI" button opens a new window for Region of Interest (ROI) analysis on the SW Group Velocity Map.

 14.   For this example, the background has a mean SW group velocity of 2.6 ± 0.3 m/s, while the inclusion presents 1.4 ± 0.2 m/s.



# How to Use (Your Own File)

To use your own data, generate a compatible .mat file containing ONLY a 3D matrix of RF data structured as (Z, X, Frames), where:

Z is the axial dimension.

X is the lateral dimension.

Frames is the temporal dimension.

Note: If the number of frames is less than 3, the file will open as RF data, but displacement calculation will fail.

The US parameters must match your file's settings: Sampling Frequency, dx, Central Frequency, Speed of Sound, and Frame Rate.

To change the default US parameter values:

Type designer in the terminal.

Open the ./gui/pyvma_ppg_white.ui file.

Edit the values as needed.
