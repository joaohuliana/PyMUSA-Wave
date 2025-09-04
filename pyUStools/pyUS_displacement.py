import math
import numpy as np
from scipy.signal import hilbert2
from scipy.signal.windows import hann

from PyQt6.QtWidgets import QApplication

from pyUStools.ui_utilities import ProgressBarWindow

import time

def framesLoupas(IQbuffer,trkprm):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Loupas's Algorithm")
    progress_bar_window.show()

    [Z, X, frames] = (IQbuffer.shape)

    c = trkprm['c']
    fc = trkprm['fc']
    winsize = trkprm['winsize']
    overlap = trkprm['overlap']
    velo = trkprm['velo']

    winoverlap = float(overlap) / 100
    winshiftpx = math.floor(winsize * (1 - winoverlap))
    if winshiftpx<=1:
        winshiftpx = 1


    Z2 = math.floor((float(Z)-(float(winsize)*winoverlap))/(float(winsize)*(1-winoverlap)))

    loupas = np.zeros((Z2+1, X, frames-1))

    I = np.zeros((winsize, X, 2))
    Q = np.zeros((winsize, X, 2))

    f = 0;
    ff = frames

    while f <= (ff-2):

        z0 = 0;
        z = 0;
        while z <= Z2 and z0 + winsize <= Z:

            if velo == 1:
                I[:, :, 0] = np.squeeze(np.real(IQbuffer[z0 : z0 + winsize, : , f]))
                I[:, :, 1] = np.squeeze(np.real(IQbuffer[z0 : z0 + winsize, :, f+1]))
                Q[:, :, 0] = np.squeeze(np.imag(IQbuffer[z0 : z0 + winsize, :, f]))
                Q[:, :, 1] = np.squeeze(np.imag(IQbuffer[z0 : z0 + winsize, :, f+1]))
            else:
                I[:, :, 0] = np.squeeze(np.real(IQbuffer[z0 : z0 + winsize, :, 0]))
                I[:, :, 1] = np.squeeze(np.real(IQbuffer[z0 : z0 + winsize, :, f + 1]))
                Q[:, :, 0] = np.squeeze(np.imag(IQbuffer[z0 : z0 + winsize, :, 0]))
                Q[:, :, 1] = np.squeeze(np.imag(IQbuffer[z0 : z0 + winsize, :, f + 1]))

            sum1 = np.zeros((1,X))
            sum2 = np.zeros((1, X))
            sum3 = np.zeros((1, X))
            sum4 = np.zeros((1, X))

            for i in range(0,winsize-1):
                sum1 = sum1 + Q[i,:,0] * I[i,:,1] - I[i,:,0]*Q[i,:,1]
                sum2 = sum2 + I[i, :, 0] * I[i, :, 1] + Q[i, :, 0] * Q[i, :, 1]
                sum3 = sum3 + Q[i, :, 0] * I[i + 1, :, 1] - I[i, :, 0] * Q[i + 1, :, 1]
                sum4 = sum4 + I[i, :, 0] * I[i + 1, :, 1] + Q[i, :, 0] * Q[i + 1, :, 1]

            #loupas[z, :, f] = (1000000) * ((c / (4 * math.pi * fc)) * (
            #        (np.arctan(sum1 / sum2)) / (1 + (np.arctan(sum3 / sum4)) / (2 * math.pi * fc))))

            loupas[z, :, f] = (1000000) * ((c / (4 * math.pi * fc)) * (
                        (np.arctan(sum1 / sum2)) / (1 + (np.arctan(sum3 / sum4)) / (2 * math.pi))))

            z0 = z0 + winshiftpx
            z = z + 1
        f = f + 1

        progress_bar_window.set_progress_value(int((f/(ff-1)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return loupas

def framesKasai(IQbuffer,trkprm):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Kasai's Algorithm")
    progress_bar_window.show()

    [Z, X, frames] = (IQbuffer.shape)
    velo = trkprm['velo']

    kasai = np.zeros((Z, X, frames-1))
    IQ1 = np.zeros((Z, X))
    IQ2 = np.zeros((Z, X))
    ImMean = np.zeros((Z, X))
    ReMean = np.zeros((Z, X))
    Im = np.zeros((Z, X))
    Re = np.zeros((Z, X))

    f = 0
    ff = frames
    while f <= (ff-2):

        if velo == 1:
            IQ1 = IQbuffer[:, :, f]
            IQ2 = IQbuffer[:, :, f + 1]
        else:
            IQ1 = IQbuffer[:, :, 0]
            IQ2 = IQbuffer[:, :, f + 1]

        ImMean = (np.imag(IQ1) + np.imag(IQ2))/2
        ReMean = (np.real(IQ1) + np.real(IQ2)) / 2

        Im = ((np.imag(IQ1) - ImMean) * (np.real(IQ2) - ReMean) -
                (np.imag(IQ2) - ImMean) * (np.real(IQ1) - ReMean))
        Re = ((np.imag(IQ1) - ImMean) * (np.imag(IQ2) - ImMean) +
              (np.real(IQ2) - ReMean) * (np.real(IQ1) - ReMean))

        kasai[:,:,f] = np.power((Im*Im + Re*Re),0.125)

        progress_bar_window.set_progress_value(int((f/(ff-1)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
        f = f + 1

    progress_bar_window.close_window()
    return kasai

def framesCrossCor(rfdata,trkprm):
    app = QApplication.instance()
    progress_bar_window = ProgressBarWindow()
    progress_bar_window.center_on_screen()
    progress_bar_window.set_title("Cross Correlation Algorithm")
    progress_bar_window.show()

    [Z, X, frames] = (rfdata.shape)

    winsize = trkprm['winsize']
    overlap = trkprm['overlap']
    velo = trkprm['velo']

    winoverlap = float(overlap) / 100
    winshiftpx = math.floor(winsize * (1 - winoverlap))

    Z2 = math.floor((float(Z)-(float(winsize)*winoverlap))/(float(winsize)*(1-winoverlap)))

    crosscor = np.zeros((Z2+1, X, frames-1))

    rf0 = np.zeros((winsize, X))
    rf1 = np.zeros((winsize, X))
    win1D = hann(winsize)
    win2D = np.tile(win1D[:, np.newaxis], (1, X))

    m = 0;
    f = 0;
    ff = frames
    noiseLevel = framesACor(rfdata[:,:,0], trkprm)

    while f <= (ff-2):

        z0 = 0;
        z = 0;
        while z <= Z2 and z0 + winsize <= Z:

            if velo == 1:
                rf0[:, :] = np.squeeze(win2D*rfdata[z0 : z0 + winsize, : , f])
                rf1[:, :] = np.squeeze(win2D*rfdata[z0 : z0 + winsize, :, f+1])
            else:
                rf0[:, :] = np.squeeze(win2D*rfdata[z0 : z0 + winsize, : , 0])
                rf1[:, :] = np.squeeze(win2D*rfdata[z0 : z0 + winsize, :, f+1])

            frf0 = np.fft.fft(rf0, axis = 0, n = 128)
            frf1 = np.fft.fft(rf1, axis = 0, n = 128)

            Rxx = np.real(np.fft.ifft(np.conj(frf0)*frf1, axis = 0));
            p = np.argmax(Rxx, axis = 0)

            for i in range(0,X):
                crosscor[z, i, f] = peakAdjust(Rxx[:,i],p[i],128-1) - noiseLevel

            z0 = z0 + winshiftpx
            z = z + 1
        f = f + 1

        progress_bar_window.set_progress_value(int((f/(ff-1)) * 100))
        progress_bar_window.repaint()
        app.processEvents()
    progress_bar_window.close_window()
    return crosscor

def framesACor(rfdata,trkprm):
    [Z, X] = (rfdata.shape)

    c = trkprm['c']
    fc = trkprm['fc']
    winsize = trkprm['winsize']
    overlap = trkprm['overlap']
    velo = trkprm['velo']

    winoverlap = float(overlap) / 100
    winshiftpx = math.floor(winsize * (1 - winoverlap))

    Z2 = math.floor((float(Z)-(float(winsize)*winoverlap))/(float(winsize)*(1-winoverlap)))

    Acor = np.zeros((Z2+1, X))

    rf0 = np.zeros((winsize, X))

    z0 = 0;
    z = 0;
    while z <= Z2 and z0 + winsize <= Z:

        rf0[:, :] = np.squeeze(rfdata[z0 : z0 + winsize, :])
        frf0 = np.fft.fft(rf0,axis=0, n=128)

        Rxx = np.real(np.fft.ifft(np.conj(frf0)*frf0,axis=0));
        p = np.argmax(Rxx, axis=0)

        for i in range(0,X):
            Acor[z, i] = peakAdjust(Rxx[:,i],p[i],128-1)

        z0 = z0 + winshiftpx
        z = z + 1

    noiseLevel = np.mean(Acor)
    return noiseLevel

def peakAdjust(Rxx, x0, N):
    perc = np.array([0, 0, 0])

    if (x0 == 0):
        xm = N
        y1 = Rxx[xm]
        xp = x0 + 1
        y3 = Rxx[xp]
        perc[0] = perc[0] + 1
        perc[2] = perc[2] + 1

    elif (x0 == N):
        xm = x0 - 1
        y1 = Rxx[xm]
        xp = 1
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

    if (x0 > N / 2):
        x0 = x0 - N

    peak = x0 + delta - 1
    return peak

