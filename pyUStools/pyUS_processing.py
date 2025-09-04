import math
import statistics, warnings
import numpy as np

from scipy.io import savemat, whosmat
from scipy.signal import hilbert2, butter, lfilter, detrend, medfilt, tukey
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.ndimage.filters import uniform_filter
from scipy.ndimage.measurements import variance
from scipy.ndimage import sobel

from PyQt6.QtWidgets import QApplication
from tkinter import filedialog as fd

from pyUStools.ui_utilities import *


#------------------------------------RF Data--------------------------------------------------
def bmode(rfdata, currentframe, compression):
    #raw = np.real(rfdata[:, :, currentframe - 1])
    raw = rfdata[:, :, currentframe - 1]
    if compression == 'Square Root (a.u.)':
        b_mode = np.sqrt(abs(hilbert2(raw)))
    elif compression == 'Logarithm (dB)':
        aux = abs(hilbert2(raw))
        aux = aux/np.amax(aux)
        b_mode = 20*np.log10(aux)
    elif compression == 'Mod Logarithm (dB)':
        aux = abs(hilbert2(raw))
        aux = aux / np.amax(aux)
        b_mode = 20 * np.log10(0.01+aux)
    else:
        b_mode = abs(hilbert2(raw))
    return b_mode

def cropRfdata(rfdata,prm,roi):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Croping RF Data")
    progress_bar_window.show()

    fs = prm['fs']
    dt = 1 / fs
    fc = prm['fc']
    df = prm['df']
    ini = prm['ini']
    end = prm['end']

    xi = roi['xi']
    xf = roi['xf']
    zi = roi['zi']
    zf = roi['zf']

    crop_rf = np.zeros((zf-zi, xf-xi, math.floor((end-ini)/df)))
    [Z, X, cframes] = (crop_rf.shape)
    [_, _, frames] = (rfdata.shape)

    i = ini
    j = 0
    while i<=(frames-1) and j<= (cframes-1):
        crop_rf[:,:,j] = rfdata[zi:zf, xi:xf, i]
        i = i + df
        j = j + 1
        progress_bar_window.set_progress_value(int((i / (frames)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()

    return crop_rf

def calculateIQ(rfdata,prm,roi):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Calculating IQ Data")
    progress_bar_window.show()

    fs = prm['fs']
    dt = 1 / fs
    fc = prm['fc']
    df = prm['df']
    ini = prm['ini']
    end = prm['end']

    xi = roi['xi']
    xf = roi['xf']
    zi = roi['zi']
    zf = roi['zf']

    crop_rf = np.zeros((zf-zi, xf-xi, math.floor((end-ini)/df)))
    [Z, X, cframes] = (crop_rf.shape)
    [_, _, frames] = (rfdata.shape)

    i = ini
    j = 0
    while i<=(frames-1) and j<= (cframes-1):
        crop_rf[:,:,j] = rfdata[zi:zf, xi:xf, i]
        i = i + df
        j = j + 1

    t = np.arange(0, Z, 1)
    t = t*dt
    t2D = np.tile(t, (X, 1)).T
    IQbuffer = np.zeros((Z,X,cframes), dtype=complex)
    for i in range(0,cframes):
        IQbuffer[:,:,i] = hilbert2(crop_rf[:,:,i])*(np.exp(-2j * np.pi * fc * t2D) / np.sqrt(2))
        progress_bar_window.set_progress_value(int((i / (cframes)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()

    return IQbuffer

def USiniavg (rfdata,N):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Average of initial frames")
    progress_bar_window.show()
    [Z, X, frames] = (rfdata.shape)
    avgUS = np.zeros((Z, X, round(frames-(N-1))))
    aux = np.zeros((Z, X))
    for i in range(0, N):
        aux = aux + rfdata[:, :, i]
    avgUS[:, :, 0] = aux / N
    count = 1
    for i in range(N, frames):
        avgUS[:, :, count] = rfdata[:, :, i]
        count = count + 1
        progress_bar_window.set_progress_value(int((i/ frames) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return avgUS

def USavg(rfdata,N):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Average of frames")
    progress_bar_window.show()
    [Z, X, frames] = (rfdata.shape)
    avgUS = np.zeros((Z, X, round(frames / N)))
    i = 0
    count = 0
    while i <= (frames-N):
        aux = np.zeros((Z, X))
        j = i
        while j <= (i + (N-1)):
            aux = aux + rfdata[:,:,j]
            j = j + 1
        avgUS[:,:,count] = aux/N
        i = i + N
        count = count + 1;
        progress_bar_window.set_progress_value(int((i / (frames-N)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return avgUS

def USinterp (rfdata):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Interpolating Channels")
    progress_bar_window.show()
    [Z, X, frames] = (rfdata.shape)
    interpUS = np.zeros((Z, 2*X, frames))
    new_x = np.linspace(0, rfdata.shape[0] - 1, 1 * rfdata.shape[0])  # Double the size along the first dimension (rows)
    new_z = np.linspace(0, rfdata.shape[1] - 1, 2 * rfdata.shape[1])  # Double the size along the second dimension (columns)

    for i in range(0, frames):
        aux = rfdata[:,:,i]
        x_values = np.arange(aux.shape[0])  # Coordinates along the first dimension
        z_values = np.arange(aux.shape[1])  # Coordinates along the second dimension

        # Create a RegularGridInterpolator
        interpolator = RegularGridInterpolator((x_values, z_values), aux, method='linear', bounds_error=False,
                                       fill_value=None)

        # Create the meshgrid for evaluation
        xx, zz = np.meshgrid(new_x, new_z, indexing='ij')

        interpUS[:,:,i] = interpolator((xx, zz))
        progress_bar_window.set_progress_value(int((i/ frames) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return interpUS

#------------------------------------Displacement Map--------------------------------------------------
def saveDispChoice():
    fdialog = SaveDisp_Dialog()
    # Define variables to store ini and end values
    choice = None

    # Connect the frames_index signal to update ini and end variables
    def update_choice(c):
        nonlocal choice
        choice = c

    fdialog.saveas.connect(update_choice)

    # Execute the dialog
    fdialog.exec()
    return choice

def saveVeloChoice():
    fdialog = SaveVelo_Dialog()
    # Define variables to store ini and end values
    choice = None

    # Connect the frames_index signal to update ini and end variables
    def update_choice(c):
        nonlocal choice
        choice = c

    fdialog.saveas.connect(update_choice)

    # Execute the dialog
    fdialog.exec()
    return choice


        #----------------------------Time Filtering------------------------------------------

def integrateDisp(dispdata):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Integrating Velocity")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    integrateDisp = np.zeros((Z, X, frames))
    for f in range(0,frames):
        integrateDisp[:,:,f] = np.sum(dispdata[:,:,0:f],axis=2)
        progress_bar_window.set_progress_value(int((f / frames) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return integrateDisp

def diffDisp(dispdata):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Diff Displacement")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    difDisp = np.zeros((Z, X, frames-1))
    for f in range(0,frames-1):
        difDisp[:,:,f] = (dispdata[:,:,f+1]-dispdata[:,:,f])/2
        progress_bar_window.set_progress_value(int((f / (frames-1)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return difDisp

def DispLPfilter(dispdata, cutoff_freq):
    if(cutoff_freq>=1):
        cutoff_freq = 0.99
    elif(cutoff_freq<=0):
        cutoff_freq = 0.01

    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Low-Pass Filter")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    pad = 50
    filtered_data = np.zeros((Z, X, frames + pad))
    aux = np.zeros((frames+pad))
    win = tukey(frames+pad, 0.5)
    # Create a low-pass Butterworth filter
    b, a = butter(3, cutoff_freq, btype='low', analog=False)
    for x in range(0, X):
        for z in range(0, Z):
            aux[int(pad/2):-int(pad/2)] = np.squeeze(dispdata[z, x, :])
            aux = aux*win
            filtered_data[z, x, :] = lfilter(b, a, aux)
        progress_bar_window.set_progress_value(int((x / (X)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    filtered_data = filtered_data[:, :, int(pad/2):-int(pad/2)]

    return filtered_data

def DispBPfilter(dispdata, cutoff_freq1, cutoff_freq2):
    if(cutoff_freq2>=1):
        cutoff_freq2 = 0.99
    elif(cutoff_freq1<=0):
        cutoff_freq1 = 0.01

    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Low-Pass Filter")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    pad = 50
    filtered_data = np.zeros((Z, X, frames + pad))
    aux = np.zeros((frames+pad))
    win = tukey(frames+pad, 0.5)
    # Create a low-pass Butterworth filter
    b, a = butter(3, [cutoff_freq1,cutoff_freq2] , btype='band', analog=False)
    for x in range(0, X):
        for z in range(0, Z):
            aux[int(pad/2):-int(pad/2)] = np.squeeze(dispdata[z, x, :])
            aux = aux*win
            filtered_data[z, x, :] = lfilter(b, a, aux)
        progress_bar_window.set_progress_value(int((x / (X)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    filtered_data = filtered_data[:, :, int(pad/2):-int(pad/2)]

    return filtered_data

def Dispinterp (dispdata):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Interpolating Frames")
    progress_bar_window.show()
    [Z, X, frames] = (dispdata.shape)
    interpDisp = np.zeros((Z, X, 2*frames))
    frames_new = np.linspace(0, frames, 2*frames)
    frames_old = np.linspace(0, frames, frames)
    for x in range(0, X):
        for z in range(0,Z):
            aux = np.squeeze(dispdata[z,x,:])
            # Create interpolation function
            interp_func = interp1d(frames_old, aux, kind='linear')
            interpaux = interp_func(frames_new)
            interpDisp[z,x,:] = interpaux
        progress_bar_window.set_progress_value(int((x/ X) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return interpDisp

def cropFrames(dispdata):
    [Z, X, frames] = (dispdata.shape)
    fdialog = Frames_Dialog()

    fdialog.sbEnd.setMaximum(frames)
    fdialog.sbEnd.setValue(frames)

    fdialog.sbIni.setMinimum(1)
    fdialog.sbIni.setValue(1)

    # Define variables to store ini and end values
    ini = None
    end = None
    # Connect the frames_index signal to update ini and end variables
    def update_ini_end(i, e):
        nonlocal ini, end
        ini, end = (i-1, e-1)

    fdialog.frames_index.connect(update_ini_end)

    # Execute the dialog
    fdialog.exec()

    if ini == -2 and end == -2:
        return dispdata
    else:
        crop = np.zeros((Z,X,(end-ini)))
        count = 0
        for f in range(ini,end):
            crop[:,:,count] = dispdata[:,:,f]
            count = count + 1
        return crop

def Dispdetrend(dispdata):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Detrending")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    filtered_data = np.zeros((Z, X, frames))
    aux_data = np.zeros((frames))
    for x in range(0, X):
        for z in range(0, Z):
            aux_data = np.squeeze(dispdata[z, x, :])
            filtered_data[z, x, :] = detrend(aux_data)
        progress_bar_window.set_progress_value(int((x / (X)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return filtered_data

        # ----------------------------Spatial Filtering------------------------------------------

def medianFilter(dispdata, mx, mz):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Median Filter")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    filtered_data = np.zeros((Z, X, frames))
    kernel_size = (mz, mx)
    for f in range(0, frames):
        filtered_data[:,:,f] = medfilt(dispdata[:,:,f], kernel_size=kernel_size)
        progress_bar_window.set_progress_value(int((f / (frames)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return filtered_data

def leeFilter(dispdata, mx, mz):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Lee's Filter")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    filtered_data = np.zeros((Z, X, frames))
    hmx = math.floor(mx/2)
    hmz = math.floor(mz/2)

    for f in range(0,frames):
        img_mean = uniform_filter(dispdata[:,:,f], (mz, mx))
        img_sqr_mean = uniform_filter(dispdata[:,:,f] ** 2, (mz, mx))
        img_variance = img_sqr_mean - img_mean ** 2

        overall_variance = variance(dispdata[:,:,f])

        #img_weights = img_variance / (img_variance + overall_variance)
        img_weights = img_variance / (img_variance + 0.5)
        filtered_data[:,:,f] = img_mean + img_weights * (dispdata[:,:,f] - img_mean)
        progress_bar_window.set_progress_value(int((f / (frames)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return filtered_data

def removeDC(dispdata):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Remove DC")
    progress_bar_window.show()

    [Z, X, frames] = (dispdata.shape)
    filtered_data = np.zeros((Z, X, frames))
    for f in range(0, frames):
        filtered_data[:, :, f] = dispdata[:, :, f] - np.mean(dispdata[:, :, f])
        progress_bar_window.set_progress_value(int((f / (frames)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return filtered_data

def directionalFilter(dispdata,dx,dz,fcoeffs,order,power,angle):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Pengfei Song Directional Filter")
    progress_bar_window.show()

    angle = np.deg2rad(angle)
    [Z, X, frames] = (dispdata.shape)

    shiftx = 1 - (X%2)/2
    shiftz = 1 - (Z%2)/2

    #Tukey window on time
    win = tukey(frames, 0.1)
    for z in range(0,Z):
        for x in range(0,X):
            aux = np.squeeze(dispdata[z,x,:])
            dispdata[z,x,:] = aux*win

    #Create ndgrid for X and Z
    lx = np.arange(X)
    lz = np.arange(Z)
    jj, ii = np.meshgrid(lx, lz)

    #Create and center grid for cossine and sine
    rhoVector = np.zeros((Z,X,2))
    rhoVector[:, :, 0] = (1/dz) * (1/Z) * (ii - int(Z/2 + shiftz))
    rhoVector[:, :, 1] = (1/dx) * (1/X) * (jj - int(X/2 + shiftx))
    rho = np.sqrt(rhoVector[:, :, 0] ** 2 + rhoVector[:, :, 1] ** 2)
    zero_indices = rho == 0
    rho[zero_indices] = 1

    #high cut and low cut coefficients in wavenumber
    hc = fcoeffs[1] * ((1/dx)/2)
    lc = fcoeffs[0] * ((1/dx)/2)

    #Create bandpass filter
    bpFilter = (1+(((hc*lc)-rho**2)/((hc-lc)*rho))**(2*order))**(-0.5)
    ind = np.isnan(bpFilter)
    bpFilter[ind] = 0

    #Calculate the filters
    dfVector = np.array([-math.sin(angle),math.cos(angle)])
    dfVector2 = np.array([math.sin(angle), -math.cos(angle)])
    normrhoVector = np.zeros((Z,X,2))
    for i in range(0,2):
        normrhoVector[:,:,i] = rhoVector[:,:,i]/rho
    ind = np.isnan(normrhoVector)
    normrhoVector[ind] = 0

    #Positive temporal frequencies
    temp = dfVector[0]*normrhoVector[:,:,0] + dfVector[1]*normrhoVector[:,:,1]
    ind = np.isnan(temp)
    temp[ind] = 0
    zero_indices = temp < 0
    temp[zero_indices] = 0

    df = (np.power(temp,power)) * bpFilter

    #Negative temporal frequencies
    temp = dfVector2[0] * normrhoVector[:, :, 0] + dfVector2[1] * normrhoVector[:, :, 1]
    ind = np.isnan(temp)
    temp[ind] = 0
    zero_indices = temp < 0
    temp[zero_indices] = 0

    dfc = (np.power(temp, power)) * bpFilter

    #Apply filter
    xf = np.fft.fft(dispdata, axis=2);
    kspdf = np.zeros((Z,X,frames),dtype=complex)

    if frames % 2 != 0:
        checkpoint = (frames + 1) / 2;
    else:
        checkpoint = frames / 2;

    #multiply the bandpass filter in kspace
    for n in range (0,frames):
        if n <= checkpoint:
            fft_1 = np.fft.fft(xf[:, :, n], axis=0)
            fft_2 = np.fft.fft(fft_1, axis=1)
            fftshift_1 = np.fft.fftshift(fft_2, axes=0)
            fftshift_2 = np.fft.fftshift(fftshift_1, axes=1)
            ksp = fftshift_2
            kspdf[:,:,n] = ksp * dfc
        else:
            fft_1 = np.fft.fft(xf[:, :, n], axis=0)
            fft_2 = np.fft.fft(fft_1, axis=1)
            fftshift_1 = np.fft.fftshift(fft_2, axes=0)
            fftshift_2 = np.fft.fftshift(fftshift_1, axes=1)
            ksp = fftshift_2
            kspdf[:, :, n] = ksp * df

        progress_bar_window.set_progress_value(int((n / (frames)) * 100))
        progress_bar_window.repaint()
        app.processEvents()

    #Inverse ffts for return from k-space
    ifftshift_1 = np.fft.ifftshift(kspdf, axes=0)
    ifftshift_2 = np.fft.ifftshift(ifftshift_1, axes=1)
    ifft_1 = np.fft.ifft(ifftshift_2, axis=0)
    ifft_2 = np.fft.ifft(ifft_1, axis=1)
    ifft_3 = np.fft.ifft(ifft_2, axis=2)
    dfiltered = np.real(ifft_3)
    progress_bar_window.close_window()

    return dfiltered

#------------------------------------Time analisys--------------------------------------------------
def powerSpectrum(signal,fs):
    # Compute the FFT of the signal
    interfft = np.fft.fft(signal, n=1024)

    # Compute the power spectrum (magnitude squared)
    interps = np.abs(interfft) ** 2

    # Frequencies corresponding to the FFT result
    # For real input, the FFT result is symmetric, and the positive frequencies are given by:
    frequencies = np.fft.fftfreq(len(interfft), d=1 / fs)
    positive_frequencies = frequencies[:len(frequencies) // 2]
    ps = interps[:len(interps) // 2]
    return ps, positive_frequencies

def areaDisp(dispmap,prm):
    app = QApplication.instance()
    area = Disp_area(dispmap,prm)
    area.exec()
    area.close()

#------------------------------------Group velocity--------------------------------------------------

def swgvelo2(data,pos_x,pos_z,dt,tlim,type):
    [_,Zm,Xm] = data.shape
    hXm = math.floor(Xm / 2);
    hZm = math.floor(Zm / 2);
    sw_pos = np.zeros((Zm,Xm))

    if type == 'max':
        A = np.squeeze(data[:,hZm,hXm])
        p0 = np.argmax(A)
        p0 = polypeak(A, p0)
        for i in range(0,Zm):
            for j in range(0,Xm):
                B = np.squeeze(data[:,i,j]);
                p1 = np.argmax(B)
                p1 = polypeak(B, p1)
                sw_pos[i,j] = (p1 - p0)

    if type == 'min':
        A = np.squeeze(data[:,hZm-1,hXm-1])
        p0 = np.argmin(A)
        p0 = polypeak(np.abs(A), p0)
        for i in range(0,Zm):
            for j in range(0,Xm):
                B = np.squeeze(data[:,i,j]);
                p1 = np.argmin(B)
                p1 = polypeak(np.abs(B), p1)
                sw_pos[i,j] = (p1 - p0)

    if type == 'corr':
        A = np.squeeze(data[:,hZm,hXm])
        A = np.fft.fft(A, axis=0, n=128)
        for i in range(0,Zm):
            for j in range(0,Xm):
                B = np.squeeze(data[:,i,j]);
                B = np.fft.fft(B, axis=0, n=128)
                C = np.real(np.fft.ifft(np.conj(A)*B, axis = 0))
                p1 = np.argmax(C)
                sw_pos[i, j] = peakAdjust(C, p1)
    gradient_x = sobel(sw_pos, axis=0)
    gradient_y = sobel(sw_pos, axis=1)

    # Compute the gradient direction in radians
    Gdir = np.arctan2(gradient_y, gradient_x)
    #dir = np.deg2rad(Gdir[hZm,hXm])
    dir = Gdir[hZm, hXm]-(math.pi/2)

    i = 0
    x_pos = np.array([])
    z_pos = np.array([])
    for k in range(-Xm,Xm):
        aux_x = round(hXm+k*math.cos(-dir))
        if(aux_x>=0 and aux_x<Xm):
            aux_z = round(hXm + k * math.sin(-dir))
            if (aux_z >= 0 and aux_z < Zm):
                x_pos = np.concatenate((x_pos,[aux_x]))
                z_pos = np.concatenate((z_pos,[aux_z]))
                i = i + 1

    t = np.array([])
    sx = np.array([])
    sz = np.array([])
    for i in range(0,len(x_pos)):
        z_aux = float(z_pos[i] - hZm)
        x_aux = float(x_pos[i] - hXm)
        nlim = math.sqrt(z_aux*z_aux+x_aux*x_aux)*tlim
        if z_pos[i] == hZm and x_pos[i] == hXm:
            t = np.concatenate((t,[sw_pos[int(z_pos[i]),int(x_pos[i])]*dt]))
            sx = np.concatenate((sx, [pos_x[int(z_pos[i]), int(x_pos[i])]]))
            sz = np.concatenate((sz, [pos_z[int(z_pos[i]), int(x_pos[i])]]))
        elif(abs(sw_pos[int(z_pos[i]),int(x_pos[i])]) >= nlim):
            t = np.concatenate((t, [sw_pos[int(z_pos[i]), int(x_pos[i])] * dt]))
            sx = np.concatenate((sx, [pos_x[int(z_pos[i]), int(x_pos[i])]]))
            sz = np.concatenate((sz, [pos_z[int(z_pos[i]), int(x_pos[i])]]))

    if len(t)>=3:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=np.RankWarning)
            coefficients = np.polyfit(t, sx, 1)
            ax = coefficients[0]  # Coefficient of x^1 (slope)
            coefficients = np.polyfit(t, sz, 1)
            az = coefficients[0]  # Coefficient of x^1 (slope)
        vx = (ax / 1000)
        vz = (az / 1000)
    else:
        vx = 0
        vz = 0
    v = np.sqrt(np.power(vx,2)+np.power(vz,2))
    return vx, vz, v

def swgvelo1(data,pos_x,pos_z,dt,tlim,type):
    [_,Zm,Xm] = data.shape
    hXm = math.floor(Xm / 2);
    hZm = math.floor(Zm / 2);
    sw_pos = np.zeros((Zm,Xm))

    if type == 'max':
        A = np.squeeze(data[:,hZm,hXm])
        p0 = np.argmax(A)
        p0 = polypeak(A, p0)
        for i in range(0,Zm):
            for j in range(0,Xm):
                B = np.squeeze(data[:,i,j]);
                p1 = np.argmax(B)
                p1 = polypeak(B, p1)
                sw_pos[i,j] = (p1 - p0)

    if type == 'min':
        A = np.squeeze(data[:,hZm-1,hXm-1])
        p0 = np.argmin(A)
        p0 = polypeak(np.abs(A), p0)
        for i in range(0,Zm):
            for j in range(0,Xm):
                B = np.squeeze(data[:,i,j]);
                p1 = np.argmin(B)
                p1 = polypeak(np.abs(B), p1)
                sw_pos[i,j] = (p1 - p0)

    if type == 'corr':
        A = np.squeeze(data[:,hZm,hXm])
        A = np.fft.fft(A, axis=0, n=128)
        for i in range(0,Zm):
            for j in range(0,Xm):
                B = np.squeeze(data[:,i,j]);
                B = np.fft.fft(B, axis=0, n=128)
                C = np.real(np.fft.ifft(np.conj(A)*B, axis = 0))
                p1 = np.argmax(C)
                sw_pos[i,j] = peakAdjust(C, p1)

    gradient_x = sobel(sw_pos, axis=0)
    gradient_y = sobel(sw_pos, axis=1)

    # Compute the gradient direction in radians
    Gdir = np.arctan2(gradient_y, gradient_x)
    #dir = np.deg2rad(Gdir[hZm,hXm])
    dir = Gdir[hZm, hXm]-(math.pi/2)

    i = 0
    x_pos = np.array([])
    z_pos = np.array([])
    for k in range(-Xm,Xm):
        aux_x = round(hXm+k*math.cos(-dir))
        if(aux_x>=0 and aux_x<Xm):
            aux_z = round(hXm + k * math.sin(-dir))
            if (aux_z >= 0 and aux_z < Zm):
                x_pos = np.concatenate((x_pos,[aux_x]))
                z_pos = np.concatenate((z_pos,[aux_z]))
                i = i + 1

    t = np.array([])
    sx = np.array([])
    sz = np.array([])
    for i in range(0, len(x_pos)):
        z_aux = float(z_pos[i] - hZm)
        x_aux = float(x_pos[i] - hXm)
        nlim = math.sqrt(z_aux * z_aux + x_aux * x_aux) * tlim
        if z_pos[i] == hZm and x_pos[i] == hXm:
            t = np.concatenate((t, [sw_pos[int(z_pos[i]), int(x_pos[i])] * dt]))
            sx = np.concatenate((sx, [pos_x[int(z_pos[i]), int(x_pos[i])]]))
            sz = np.concatenate((sz, [pos_z[int(z_pos[i]), int(x_pos[i])]]))
        elif (abs(sw_pos[int(z_pos[i]), int(x_pos[i])]) >= nlim):
            t = np.concatenate((t, [sw_pos[int(z_pos[i]), int(x_pos[i])] * dt]))
            sx = np.concatenate((sx, [pos_x[int(z_pos[i]), int(x_pos[i])]]))
            sz = np.concatenate((sz, [pos_z[int(z_pos[i]), int(x_pos[i])]]))
    return t, sx, sz

def areaVelo(velomap,prm):
    app = QApplication.instance()
    area = Velo_area(velomap,prm)
    area.exec()
    area.close()


#------------------------------------Phase velocity--------------------------------------------------
def phaseVelo(dispmap,prm):
    app = QApplication.instance()
    phase = phase_velo(dispmap,prm)
    dispdata = None
    freqdata = None

    # Connect the frames_index signal to update ini and end variables
    def update_curve(f,d):
        nonlocal dispdata, freqdata
        dispdata = d
        freqdata = f

    phase.dispcurvedata.connect(update_curve)

    phase.exec()
    phase.close()
    return freqdata,dispdata


#------------------------------------Utilities--------------------------------------------------
def mirror_border(matrix, m_x, m_y):
    rows, cols = matrix.shape
    mirrored_matrix = np.zeros((rows + 2 * m_x, cols + 2 * m_y), dtype=matrix.dtype)
    mirrored_matrix[m_x:rows + m_x, m_y:cols + m_y] = matrix

    # Top border
    for i in range(m_x):
        mirrored_matrix[i, :] = mirrored_matrix[2 * m_x - i - 1, :]

    # Bottom border
    for i in range(m_x):
        mirrored_matrix[rows + m_x + i, :] = mirrored_matrix[rows + m_x - i - 1, :]

    # Left border
    for i in range(m_y):
        mirrored_matrix[:, i] = mirrored_matrix[:, 2 * m_y - i - 1]

    # Right border
    for i in range(m_y):
        mirrored_matrix[:, cols + m_y + i] = mirrored_matrix[:, cols + m_y - i - 1]

    return mirrored_matrix

def polypeak(data,x1):
    N = len(data);
    if (x1 == 0):
        peak = x1
        return peak
    elif(x1 == N-1):
        peak = x1
        return peak
    else:
        x0 = x1 - 1
        x2 = x1 + 1
        a = 2 * (2 * data[x1] - data[x2] - data[x0])
        b = data[x2] - data[x0]

        if (a == 0):
            delta = 0
        else:
            delta = b / a

        peak = x1 + delta
        return peak


def peakAdjust(Rxx, x0):
    N = len(Rxx);
    perc = np.array([0, 0, 0])

    if (x0 == 0):
        xm = N - 1
        y1 = Rxx[xm]
        xp = x0 + 1
        y3 = Rxx[xp]
        perc[0] = perc[0] + 1
        perc[2] = perc[2] + 1

    elif (x0 == N-1):
        xm = x0 - 1
        y1 = Rxx[xm]
        xp = 0
        y3 = Rxx[xp]
        perc[1] = perc[1] + 1
        perc[2] = perc[2] + 1

    else:
        xm = x0 - 1
        xp = x0 + 1
        y1 = Rxx[xm]
        y3 = Rxx[xp]
        perc[2] = perc[2] + 1

    a = 2 * (2 * Rxx[x0] - Rxx[xp] - Rxx[xm])
    b = Rxx[xp] - Rxx[xm]

    if (a == 0):
        delta = 0
    else:
        delta = b / a

    if (x0 > (N+1) / 2):
        x0 = x0 - N

    #peak = x0 + delta - 1
    peak = x0 + delta
    return peak



def fill_outliers(data, method='median', threshold=1.5):
    # Make a copy of the original data
    filled_data = data.copy()

    # Calculate the lower and upper quantiles
    q1, q3 = np.percentile(data, [25, 75])

    # Calculate the IQR (Interquartile Range)
    iqr = q3 - q1

    # Define the lower and upper bounds
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr

    # Replace outliers based on the method
    if method == 'median':
        outlier_values = np.where((data < lower_bound) | (data > upper_bound))[0]
        filled_data[outlier_values] = np.median(data)
    elif method == 'mean':
        outlier_values = np.where((data < lower_bound) | (data > upper_bound))[0]
        filled_data[outlier_values] = np.mean(data)
    else:
        raise ValueError("Invalid method. Choose 'median' or 'mean'.")

    return filled_data