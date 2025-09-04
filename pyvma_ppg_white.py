import sys, time
import os
import math
import ctypes
from pathlib import Path

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6 import uic
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtWidgets import QVBoxLayout, QHBoxLayout,QApplication, QMainWindow, QSizePolicy
from PyQt6.QtGui import QScreen, QFont
from tkinter import filedialog as fd
import numpy as np
import scipy
from scipy.io import savemat, whosmat
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import Rectangle, RectangleSelector
import matplotlib.ticker as ticker

#pyUStools
from pyUStools.pyUS_processing import *
from pyUStools.pyUS_displacement import *
from pyUStools.ui_utilities import Tracking_Dialog, Gvelocity_Dialog, ProgressBarWindow

class pyMMUS(QtWidgets.QMainWindow):
# ---------------------------------------------General UI---------------------------------------------------------------
    def __init__(self):
        super().__init__()

        # Get the path to the bundled UI file
        ui_file_path = self.resource_path("gui/pyvma_ppg_white.ui")
        # Load the UI file
        uic.loadUi(ui_file_path, self)

        # pyinstaller --onefile --add-data ".\gui\pyvma_ppg_white.ui:gui" --add-binary '.\C_functions\cmmus.dll;.' --noconsole pyvma_ppg_white.py
        # pip freeze > requirements.txt
        # pip install -r requirements.txt

        #self.setWindowFlag(Qt.WindowType.WindowMaximizeButtonHint, False)
        #self.setFixedSize(1366, 700)
        self.setMinimumSize(1366,700)
        self.setSizeIncrement(1,1)
        primary_screen = app.primaryScreen()
        screen_geometry = primary_screen.geometry()
        self.screenW = screen_geometry.width()
        self.screenH = screen_geometry.height()
        self.initialResizeEvent = False

        #center mainwindow
        qr = self.frameGeometry()
        cp = self.screen().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

        self.initUI()

    def resource_path(self, relative_path):
        # Function to get the path to a bundled resource
        if getattr(sys, 'frozen', False):
            # Running in a bundle (PyInstaller)
            return os.path.join(sys._MEIPASS, relative_path)
        else:
            # Running in a normal Python environment
            return os.path.join(os.path.abspath("."), relative_path)

    def initUI(self):
        self.menuExit.triggered.connect(self.close)
        #BMaxis
        layout = QVBoxLayout(self.bmAxis)
        layout.setContentsMargins(1,1,1,1)
        layout.addStretch()
        self.bmfig = plt.Figure()
        self.bmfig.set_facecolor((0.0, 0.0, 0.0))
        self.bmax = self.bmfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.bmax.set_facecolor((0, 0, 0))
        self.bmax.set_xlabel('Lateral [mm]', fontsize=9)
        self.bmax.set_ylabel('Axial [mm]', fontsize=9)
        self.bmax.tick_params(axis='x', labelsize=9)
        self.bmax.tick_params(axis='y', labelsize=9)
        self.bmcanvas = FigureCanvas(self.bmfig)
        layout.addWidget(self.bmcanvas)
        layout.addStretch()

        #Disp Axis
        displayout = QVBoxLayout(self.dispAxis)
        displayout.setContentsMargins(1, 1, 1, 1)
        displayout.addStretch()
        self.dispfig = plt.Figure()
        self.dispfig.set_facecolor((0.0, 0.0, 0.0))
        self.dispax = self.dispfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.dispax.set_facecolor((0, 0, 0))
        self.dispax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispax.tick_params(axis='x', labelsize=9)
        self.dispax.tick_params(axis='y', labelsize=9)
        self.dispcanvas = FigureCanvas(self.dispfig)
        displayout.addWidget(self.dispcanvas)
        displayout.addStretch()

        #Graph axis
        graphlayout = QVBoxLayout(self.graphAxis)
        graphlayout.setContentsMargins(1, 1, 1, 1)
        graphlayout.addStretch()
        self.graphfig = plt.Figure()
        self.graphfig.set_facecolor((0.0, 0.0, 0.0))
        self.graphax = self.graphfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.graphax.set_facecolor((0, 0, 0))
        self.graphax.set_xlabel('X', fontsize=9)
        self.graphax.set_ylabel('Y', fontsize=9)
        self.graphax.tick_params(axis='x', labelsize=9)
        self.graphax.tick_params(axis='y', labelsize=9)
        self.graphcanvas = FigureCanvas(self.graphfig)
        graphlayout.addWidget(self.graphcanvas)
        graphlayout.addStretch()

        self.bmloaded = False

        # parameters
        self.c = self.sbC.value()
        self.dx = self.sbdx.value()
        self.fs = (self.sbSf.value() * 1000000)
        self.fc = (self.sbFc.value() * 1000000)
        self.dz = (self.c / (2 * self.fs)) * 1000
        self.fr = self.sbFr.value()
        self.dt = float(1 / self.fr)

        self.menuLoadmat.triggered.connect(self.load_mat)
        self.menuClosescene.triggered.connect(self.close_mat)
        self.slbm.valueChanged.connect(self.bmslider_changed)
        self.pbLoadprm.clicked.connect(self.load_prm)
        self.pbSaveBM.clicked.connect(self.saveBMpng)
        self.pbSaveBM.setVisible(False)
        self.pbSaveDisp.clicked.connect(self.saveDispmat)
        self.pbSaveDisp.setVisible(False)
        self.cbCompression.activated.connect(self.compression_type)

        #ROI Panel
        self.pbSROI.clicked.connect(self.select_roi)
        self.pbResetROI.clicked.connect(self.reset_roi)
        self.cbROItype.activated.connect(self.roi_type)

        self.sbXi.valueChanged.connect(self.roi_pbchanged)
        self.sbXf.valueChanged.connect(self.roi_pbchanged)
        self.sbZi.valueChanged.connect(self.roi_pbchanged)
        self.sbZf.valueChanged.connect(self.roi_pbchanged)

        #US Processing Panel
        self.pbUSavg.clicked.connect(self.US_avg)
        self.pbUSiniavg.clicked.connect(self.US_iniavg)
        self.pbUSinterp.clicked.connect(self.US_interp)

        #US Parameter Panel
        self.sbSf.valueChanged.connect(self.sf_changed)
        self.sbdx.valueChanged.connect(self.dx_changed)
        self.sbC.valueChanged.connect(self.c_changed)
        self.sbFr.valueChanged.connect(self.fr_changed)
        self.sbFc.valueChanged.connect(self.fc_changed)

        #Displacement
        self.pbLoupas.clicked.connect(self.loupas)
        self.pbKasai.clicked.connect(self.kasai)
        self.pbCrossCor.clicked.connect(self.cross_cor)
        self.pbTrkset.clicked.connect(self.tracking_prm)
        self.slDisp.valueChanged.connect(self.dispslider_changed)
        self.slDispMax.valueChanged.connect(self.dispslider_changed)
        self.slDispMax.setVisible(False)
        self.lbDispMax.setVisible(False)
        self.slDispMin.valueChanged.connect(self.dispslider_changed)
        self.slDispMin.setVisible(False)
        self.lbDispMin.setVisible(False)

        #Time filtering
        self.pbLPfilter.clicked.connect(self.LPfilter)
        self.pbBPfilter.clicked.connect(self.BPfilter)
        self.pbDispinterp.clicked.connect(self.interpDisp)
        self.pbSelectFrames.clicked.connect(self.frames_selection)
        self.pbDetrend.clicked.connect(self.detrendDisp)

        #Spatial filtering
        self.pbMedianfilt.clicked.connect(self.medianDisp)
        self.pbLeefilt.clicked.connect(self.leeDisp)
        self.pbFrostfilt.clicked.connect(self.frostDisp)
        self.pbRemoveDC.clicked.connect(self.dcDisp)
        self.pbDirfilt.clicked.connect(self.dircall)
        self.pbFiltprm.clicked.connect(self.filter_prm)

        #Time Analisys
        self.pbSignal.clicked.connect(self.signal_point)
        self.pbPSD.clicked.connect(self.signal_PSD)
        self.pbArea.clicked.connect(self.areaAnalysis)

        #G-Velocity
        self.pbGVelo2D.clicked.connect(self.ToF)
        self.pbGVelo.clicked.connect(self.sToF_point)
        self.pbGVeloSet.clicked.connect(self.gvelo_set)
        self.slVeloMax.valueChanged.connect(self.veloslider_changed)
        self.slVeloMax.setVisible(False)
        self.lbVeloMax.setVisible(False)
        self.slVeloMin.valueChanged.connect(self.veloslider_changed)
        self.slVeloMin.setVisible(False)
        self.lbVeloMin.setVisible(False)
        self.cbVeloshow.setVisible(False)
        self.cbVeloshow.setVisible(False)
        self.cbVeloshow.activated.connect(self.veloshow_change)
        self.lbVeloval.setVisible(False)
        self.pbVelocursor.setVisible(False)
        self.pbVelocursor.clicked.connect(self.velo_cursor)
        self.pbSaveVelo.clicked.connect(self.saveVelomat)
        self.pbSaveVelo.setVisible(False)
        self.pbROIVelo.clicked.connect(self.veloROI)
        self.pbROIVelo.setVisible(False)

        #P-Velocity
        self.pbPVelo.clicked.connect(self.phaseVelo1)
        self.pbPVeloSet.clicked.connect(self.pvelo_prm)




        self.rect = None
        self.ROItext = None
        self.Signaltext = None
        self.Signalmarker = None
        #Tracking
        self.winsize = 20
        self.overlap = 75
        self.inif = 2
        self.skipf = 1
        #filter
        self.mx = 3
        self.mz = 3
        self.angle = 0
        self.lc = 0.01
        self.hc = 3
        self.order = 3
        self.power = 0.05
        #G-velocity
        self.dist = 1.2
        self.N = int(5)
        self.gtr = 20
        self.velocursorflag = False
        #P-velocity
        self.zeropad = 2048
        self.pvelotr = -6
        self.pvelofc = 1000
        self.pvelodf = 25
        self.pvelodensity = 1000

    def activate_guiBM(self, value):

        self.menuClosescene.setEnabled(value)

        #BMode axis
        self.bmfig.set_facecolor((0.94, 0.94, 0.94))
        self.cbCompression.setEnabled(value)
        self.slbm.setEnabled(value)
        self.pbSaveBM.setEnabled(value)
        self.pbSaveBM.setVisible(value)


        #US ROI
        self.pbSROI.setEnabled(value)
        self.sbXi.setEnabled(value)
        self.sbXf.setEnabled(value)
        self.sbZi.setEnabled(value)
        self.sbZf.setEnabled(value)
        self.cbROItype.setEnabled(value)

        #US processing
        self.pbUSavg.setEnabled(value)
        self.sbUSavg.setEnabled(value)
        self.pbUSiniavg.setEnabled(value)
        self.sbUSiniavg.setEnabled(value)
        self.pbUSinterp.setEnabled(value)

        #Displacement
        self.pbLoupas.setEnabled(value)
        self.pbKasai.setEnabled(value)
        self.pbCrossCor.setEnabled(value)
        self.cbParticle.setEnabled(value)
        self.pbTrkset.setEnabled(value)

    def activate_guiDisp(self, value):
        #Disp axis
        self.dispfig.set_facecolor((0.94, 0.94, 0.94))
        self.slDisp.setEnabled(value)
        self.slDispMax.setEnabled(value)
        self.slDispMax.setVisible(value)
        self.lbDispMax.setEnabled(value)
        self.lbDispMax.setVisible(value)

        self.slDispMin.setEnabled(value)
        self.slDispMin.setVisible(value)
        self.lbDispMin.setEnabled(value)
        self.lbDispMin.setVisible(value)

        self.pbSaveDisp.setEnabled(value)
        self.pbSaveDisp.setVisible(value)

        #Time filtering
        self.pbLPfilter.setEnabled(value)
        self.sbLPfilter.setEnabled(value)
        self.pbBPfilter.setEnabled(value)
        self.sbBPfilter1.setEnabled(value)
        self.pbDispinterp.setEnabled(value)
        self.pbSelectFrames.setEnabled(value)
        self.pbDetrend.setEnabled(value)

        #Spatial Filtering
        self.pbMedianfilt.setEnabled(value)
        self.pbLeefilt.setEnabled(value)
        self.pbFrostfilt.setEnabled(value)
        self.pbDirfilt.setEnabled(value)
        self.pbRemoveDC.setEnabled(value)
        self.pbFiltprm.setEnabled(value)

        #Time analysis
        self.pbSignal.setEnabled(value)
        self.pbArea.setEnabled(value)

        #G-Velocity
        self.pbGVelo2D.setEnabled(value)
        self.pbGVelo.setEnabled(value)
        self.pbGVeloSet.setEnabled(value)
        self.cbGmethod.setEnabled(value)

        #P-Velocity
        self.pbPVelo.setEnabled(value)
        self.pbPVeloSet.setEnabled(value)

    def activate_velosliders(self, value):
        self.slVeloMax.setVisible(value)
        self.lbVeloMax.setVisible(value)
        self.slVeloMin.setVisible(value)
        self.lbVeloMin.setVisible(value)
        self.cbVeloshow.setVisible(value)
        self.pbVelocursor.setVisible(value)

        self.slVeloMax.setEnabled(value)
        self.slVeloMin.setEnabled(value)
        self.cbVeloshow.setEnabled(value)
        self.pbVelocursor.setEnabled(value)
        self.pbSaveVelo.setEnabled(value)
        self.pbSaveVelo.setVisible(value)
        self.pbROIVelo.setEnabled(value)
        self.pbROIVelo.setVisible(value)


    def resizeEvent(self, event):

        guiObj = [self.bmAxis, self.dispAxis, self.graphAxis,
                          self.panelUSroi, self.panelUSproc, self.panelUSprm, self.panelDisp,
                          self.panelTimefilt, self.panelSpacefilt, self.panelTimeanalysis,
                          self.panelSWGroup, self.panelSWPhase, self.lbCompression, self.cbCompression,
                          self.pbSaveBM, self.slbm, self.lbFramebm, self.pbSaveDisp, self.lbFileName,
                          self.slDisp, self.lbFrameDisp, self.slDispMin, self.slDispMax, self.lbDispMax,
                          self.lbDispMin, self.pbVelocursor, self.lbVeloval, self.cbVeloshow,
                          self.slVeloMax, self.slVeloMin, self.lbVeloMax, self.lbVeloMin, self.pbROIVelo,
                          self.pbSaveVelo]
        RoiObj = [self.lbXi, self.sbXi, self.lbZi, self.sbZi, self.sbXf, self.lbXf, self.sbZf, self.lbZf,
                          self.lbROItext, self.cbROItype, self.pbSROI, self.pbResetROI]
        USprocObj = [self.pbUSavg, self.pbUSiniavg, self.sbUSavg, self.sbUSiniavg,
                             self.lbUSinterp, self.pbUSinterp]
        USprmObj = [self.lbSf, self.sbSf, self.lbSf_u, self.lbFc, self.sbFc, self.lbCf_u,
                            self.lbdx, self.sbdx, self.lbdx_u, self.lbFr, self.sbFr, self.lbFr_u,
                            self.lbC, self.sbC, self.lbC_u, self.lbFr, self.pbLoadprm]
        TrkprmObj = [self.pbCrossCor, self.pbLoupas, self.pbKasai, self.lbParticle,
                            self.cbParticle, self.pbTrkset]
        TimefiltObj = [self.pbLPfilter, self.pbBPfilter, self.sbBPfilter1, self.sbLPfilter, self.pbDispinterp,
                       self.pbSelectFrames, self.pbDetrend, self.lbLPunity, self.lbBPunity]
        SpacefiltObj = [self.pbLeefilt, self.pbFrostfilt, self.pbMedianfilt, self.pbRemoveDC,
                               self.pbDirfilt, self.pbFiltprm]
        TimeanaObj = [self.pbSignal, self.pbArea, self.pbPSD]
        GveloObj = [self.pbGVelo2D, self.pbGVelo, self.lbGmethod, self.cbGmethod,
                                self.pbGVeloSet]
        PveloObj = [self.pbPVelo, self.pbPVeloSet]

        if self.initialResizeEvent:

                w = self.size().width() / 1366
                h = self.size().height() / 700
                lw = math.sqrt(self.size().width() ** 2 + self.size().height() ** 2) / math.sqrt(1366 ** 2 + 700 ** 2)
                font = QFont()
                font.setPointSize(int(8*lw))
                for obj in guiObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in RoiObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in USprocObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in USprmObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in TrkprmObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in TimefiltObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in SpacefiltObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in TimeanaObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in GveloObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)
                for obj in PveloObj:
                    cur = obj.origeometry
                    obj.setGeometry(int(cur.x() * w), int(cur.y() * h), int(cur.width() * w), int(cur.height() * h))
                    obj.setFont(font)

                #self.bmcanvas.setGeometry(int(self.bmcanvasgeo.x()*w), int(self.bmcanvasgeo.y()*h),
                #                          int(self.bmcanvasgeo.width()*w), int(self.bmcanvasgeo.height()*h))
                self.bmcanvas.setGeometry(round(self.bmcanvasgeo.x()*w), round(self.bmcanvasgeo.y()*h),
                                              round(self.bmcanvasgeo.width()*w), round(self.bmcanvasgeo.height()*h))

                self.dispcanvas.setGeometry(int(self.bmcanvasgeo.x() * w), int(self.bmcanvasgeo.y() * h),
                                              int(self.bmcanvasgeo.width() * w), int(self.bmcanvasgeo.height() * h))
                self.graphcanvas.setGeometry(int(self.bmcanvasgeo.x() * w), int(self.bmcanvasgeo.y() * h),
                                              int(self.bmcanvasgeo.width() * w), int(self.bmcanvasgeo.height() * h))

                self.bmfig.set_size_inches([self.bmfiggeo[0] * w, self.bmfiggeo[1] * h])
                self.dispfig.set_size_inches([self.bmfiggeo[0] * w, self.bmfiggeo[1] * h])
                self.graphfig.set_size_inches([self.bmfiggeo[0] * w, self.bmfiggeo[1] * h])
                #self.dispfig.set_size_inches(a[0] * w, a[1] * h)
                #self.graphfig.set_size_inches(a[0] * w, a[1] * h)
                self.bmax.set_position = [0.17, 0.14, 0.75, 0.8]


        else:
            self.initialResizeEvent = True
            for obj in guiObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in RoiObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in USprocObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in USprmObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in TrkprmObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in TimefiltObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in SpacefiltObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in TimeanaObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in GveloObj:
                setattr(obj, 'origeometry', obj.geometry())
            for obj in PveloObj:
                setattr(obj, 'origeometry', obj.geometry())
            self.bmfiggeo = [4.1, 3]
            self.bmcanvasgeo = self.bmcanvas.geometry()

# -----------------------------------------------Load mat---------------------------------------------------------------
    def load_mat(self):
        if self.bmloaded == True:
            self.close_mat()

        filetypes = (
            ('Mat files', '*.mat'),
        )
        matfilename = fd.askopenfilename(
            title='Select a .mat file',
            initialdir='/',
            filetypes=filetypes)
        if matfilename:
            #load
            auxnamelist = whosmat(matfilename)
            auxname = auxnamelist[0]
            auxload = scipy.io.loadmat(matfilename)
            self.rfdata = auxload[auxname[0]]
            if self.rfdata.shape == (1,1):
                dialog = Error_Dialog("Load Error")
                # Set the label text
                dialog.set_text("Load error! The .mat file must contain only a 3D numeric matrix.")
                # Show the dialog
                dialog.exec()
                return
            if self.rfdata.ndim < 3:
                dialog = Error_Dialog("Frame error")
                # Set the label text
                dialog.set_text("Load error! RF data has just one frame")
                # Show the dialog
                dialog.exec()
                return

            filename_only = matfilename.split("/")[-1]
            filename_only =filename_only.replace(".mat", "")
            filename_only = filename_only.replace("_", " ")

            self.lbFileName.setText(str(filename_only))
            self.lbFileName.setAlignment(Qt.AlignmentFlag.AlignHCenter)  # Center the text horizontally
            self.lbFileName.setStyleSheet("font-weight: bold; font-size: 12pt;")

            #dimensions
            [self.Z, self.X, self.frames] = (self.rfdata.shape)
            self.endf = self.frames
            self.programmatic_change = True
            self.rfdata_frames_update()
            self.reset_roi()
            self.programmatic_change = False

            #parameters
            self.c = self.sbC.value()
            self.dx = self.sbdx.value()
            self.fs = (self.sbSf.value() * 1000000)
            self.fc = (self.sbFc.value() * 1000000)
            self.dz = (self.c / (2 * self.fs)) * 1000
            self.fr = self.sbFr.value()
            self.dt = float(1 / self.fr)

            self.bmloaded = True

            self.activate_guiBM(True)
            self.update_frame()
            self.bmfiggeo = self.bmfig.get_size_inches()
        else:
            return

    def rfdata_frames_update(self):
        self.slbm.setMinimum(1)  # Set the minimum value
        self.slbm.setMaximum(self.frames)
        self.endf = self.frames
        self.inif = 2

    def close_mat(self):
        self.activate_guiBM(False)
        self.activate_guiDisp(False)
        self.activate_velosliders(False)
        if self.velocursorflag == True:
            self.velo_cursor()
        self.bmax.clear()
        self.bmfig.set_facecolor((0, 0, 0))
        self.bmax.set_facecolor((0, 0, 0))
        self.bmcanvas.draw()
        self.dispax.clear()
        self.dispfig.set_facecolor((0, 0, 0))
        self.dispax.set_facecolor((0, 0, 0))
        self.dispcanvas.draw()
        self.graphax.clear()
        self.graphfig.set_facecolor((0, 0, 0))
        self.graphax.set_facecolor((0, 0, 0))
        self.graphcanvas.draw()
        self.bmloaded = False

    def saveBMpng(self):
        filetypes = (('PNG files', '*.png'),)

        # Open the file dialog for saving
        filename = fd.asksaveasfilename(
            title='Save As',
            filetypes=filetypes,
            defaultextension='.png'  # Optional: specify default extension
        )

        # Check if a filename was selected
        if filename:
            fig = self.bmax.get_figure()
            fig.savefig(filename)
        else:
            return

    def saveDispmat(self):
        choice = saveDispChoice()
        if choice == 'Data':
            filetypes = (
                ('Mat files', '*.mat'),
            )

            # Open the file dialog for saving
            filename = fd.asksaveasfilename(
                title='Save As',
                filetypes=filetypes,
                defaultextension='.mat'  # Optional: specify default extension
            )

            # Check if a filename was selected
            if filename:
                disp_map = self.dispmap
                dz = self.dispdz
                dx = self.dispdx
                dt = self.dt
                # Create a dictionary with key-value pairs for each variable
                variables_dict = {
                    'disp_map': disp_map,
                    'dz': dz,
                    'dx': dx,
                    'dt': dt
                }

                # Save the variables as a .mat file
                savemat(filename, variables_dict)
        elif choice == 'Figure':
            filetypes = (('PNG files', '*.png'),)

            # Open the file dialog for saving
            filename = fd.asksaveasfilename(
                title='Save As',
                filetypes=filetypes,
                defaultextension='.png'  # Optional: specify default extension
            )

            # Check if a filename was selected
            if filename:
                fig = self.dispax.get_figure()
                fig.savefig(filename,dpi=300)
            else:
                return
        else:
            return

    def saveVelomat(self):
        choice = saveVeloChoice()
        if choice == 'Data':
            filetypes = (
                ('Mat files', '*.mat'),
            )

            # Open the file dialog for saving
            filename = fd.asksaveasfilename(
                title='Save As',
                filetypes=filetypes,
                defaultextension='.mat'  # Optional: specify default extension
            )

            # Check if a filename was selected
            if filename:
                veloshow = self.cbVeloshow.currentText()
                if veloshow == 'Vx':
                    velo_map =self.vx
                elif veloshow == 'Vz':
                    velo_map = self.vz
                elif veloshow == 'Velocity':
                    velo_map = self.v

                dz = self.dispdz
                dx = self.dispdx
                dt = self.dt
                # Create a dictionary with key-value pairs for each variable
                variables_dict = {
                    'velo_map': velo_map,
                    'dz': dz,
                    'dx': dx,
                    'dt': dt
                }

                # Save the variables as a .mat file
                savemat(filename, variables_dict)
        elif choice == 'Figure':
            filetypes = (('PNG files', '*.png'),)

            # Open the file dialog for saving
            filename = fd.asksaveasfilename(
                title='Save As',
                filetypes=filetypes,
                defaultextension='.png'  # Optional: specify default extension
            )

            # Check if a filename was selected
            if filename:
                fig = self.graphax.get_figure()
                fig.savefig(filename,dpi=300)
            else:
                return
        else:
            return

    def load_prm(self):
        filetypes = (
        ('Mat files', '*.mat'),
        )
        prmfilename = fd.askopenfilename(
            title='Select a .mat file',
            initialdir='/',
            filetypes=filetypes)
        if prmfilename:
            # load
            auxnamelist = whosmat(prmfilename)
            auxname = auxnamelist[0]
            mat_data = scipy.io.loadmat(prmfilename)
            prms = mat_data[auxname[0]]
            for var_name, var_value in mat_data.items():
                if isinstance(var_value, np.ndarray):
                    # Set the flag to indicate that at least one structure was found
                    structures_found = True

                    # Access the structure
                    struct_data = var_value

                    # Get the names of the fields
                    field_names = struct_data.dtype.names
                    # Access the values corresponding to each field
                    if field_names is not None:
                        for field_name in field_names:
                            field_value = struct_data[field_name]
                            if field_name == 'fs':
                                value = field_value[0][0][0][0]
                                self.fs = (value) / 1000000
                                self.sbSf.setValue(self.fs)
                            elif field_name == 'fc':
                                value = field_value[0][0][0][0]
                                self.fc = (value) / 1000000
                                self.sbFc.setValue(self.fc)
                            elif field_name == 'c':
                                value = field_value[0][0][0][0]
                                self.c = int(value)
                                self.sbC.setValue(self.c)
                            elif field_name == 'frameRate':
                                value = field_value[0][0][0][0]
                                self.fr = (value)
                                self.sbFr.setValue(self.fr)

#------------------------------------------------BMode Image------------------------------------------------------------

    def show_image(self):
        self.bmax.clear()
        compression = self.cbCompression.currentText()
        b_mode = bmode(self.rfdata, self.current_frame, compression);

        extent = [0, self.X * self.dx, self.Z * self.dz,0]
        self.bmax.imshow(b_mode, cmap='gray',
                        extent=extent)

        self.bmax.set_xlabel('Lateral [mm]', fontsize=9)
        self.bmax.set_ylabel('Axial [mm]', fontsize=9)
        self.bmax.tick_params(axis='x', labelsize=9)
        self.bmax.tick_params(axis='y', labelsize=9)
        self.bmcanvas.draw()
        self.generate_roi()

    def update_frame(self):
        self.current_frame = self.slbm.value()
        self.lbFramebm.setText(f'{self.current_frame} / {self.frames} ')
        self.show_image()

    def bmslider_changed(self):
        if self.programmatic_change:
            return
        else:
            self.update_frame()

    def compression_type(self):
        if self.bmloaded == True:
            self.update_frame()

#-------------------------------------------------Disp Image------------------------------------------------------------
    def dispshow_image(self):
        self.dispax.clear()
        min = float(self.slDispMin.value())/100
        max = float(self.slDispMax.value())/100
        extent = [self.roixi, self.roixf, self.roizf, self.roizi]
        image = self.dispax.imshow(self.dispmap[:,:,self.dispcurrent_frame-1], cmap='jet',
                        extent=extent,vmin=min, vmax=max,interpolation='none')

        self.dispax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispax.tick_params(axis='x', labelsize=9)
        self.dispax.tick_params(axis='y', labelsize=9)
        self.dispcanvas.draw()
        self.dispax.set_aspect('equal')

    def dispupdate_frame(self):
        min = float(self.slDispMin.value()) / 100
        max = float(self.slDispMax.value()) / 100
        self.dispcurrent_frame = self.slDisp.value()
        self.lbFrameDisp.setText(f'{self.dispcurrent_frame} / {self.dispframes} ')
        self.lbDispMax.setText(f'{max: .2f}')
        self.lbDispMin.setText(f'{min: .2f}')
        self.dispshow_image()

    def dispslider_changed(self):
        self.dispupdate_frame()

# ------------------------------------------------BMode ROI-------------------------------------------------------------
    def select_roi(self):
        if self.rect is not None:
            self.rect.remove()
            self.rect = None
            self.bmcanvas.draw()
        if self.ROItext is not None:
            self.ROItext.remove()
            self.bmcanvas.draw()
        self.ROItext = self.bmax.text(0.5, 0.5, "Select a ROI \n Right-Click to finish",
                                    color='red', fontsize=12, ha='center', va='center', transform=self.bmax.transAxes, )
        self.bmcanvas.draw()
        self.rs = RectangleSelector(self.bmax, self.draw_roi, useblit=True,
                                    button=[1], minspanx=5, minspany=5, spancoords='pixels',
                                    interactive=True,  props = dict(facecolor='green', edgecolor='green', alpha=0.3, fill=True))
        self.is_drawing = True
        self.connection_id = self.bmfig.canvas.mpl_connect('button_press_event', self.finish_roi)
        self.pbSROI.setEnabled(False)
        self.activate_guiBM(False)
        self.panelUSprm.setEnabled(False)
        self.bmAxis.setCursor(Qt.CursorShape.CrossCursor)

    def draw_roi(self, eclick, erelease, *args, **kwargs):
        if self.is_drawing:
            if self.rect is not None:
                self.rect.remove()
                self.rect = None
                self.bmcanvas.draw()
            # Draw a red rectangle on the selected region
            self.rect = Rectangle((eclick.xdata, eclick.ydata),
                             erelease.xdata - eclick.xdata,
                             erelease.ydata - eclick.ydata,
                             linewidth=2, edgecolor='red', facecolor='none', linestyle = 'dashed')
            self.bmax.add_patch(self.rect)
            self.bmcanvas.draw()
            #

    def finish_roi(self,event):
        #if event.dblclick and event.button == 1:
        if event.button == 3:
            if self.rect is not None:
                self.rs.set_visible(False)
                self.rs.set_active(False)
                self.rs.disconnect_events()
                self.rs = None
                self.pbSROI.setEnabled(True)
                self.pbResetROI.setEnabled(True)
                self.roixi = self.rect.get_x()
                self.roizi = self.rect.get_y()
                self.roixf = self.roixi + self.rect.get_width()
                self.roizf = self.roizi + self.rect.get_height()
                self.ROItext.remove()
                self.ROItext = None
                self.rect.remove()
                self.rect = None
                self.bmcanvas.draw()
                self.bmfig.canvas.mpl_disconnect(self.connection_id)
                self.is_drawing = False
                self.programmatic_change = True
                self.sbXi.setValue(self.roixi)
                self.sbXf.setValue(self.roixf)
                self.sbZi.setValue(self.roizi)
                self.sbZf.setValue(self.roizf)
                self.programmatic_change = False
                self.generate_roi()
                self.activate_guiBM(True)
                self.panelUSprm.setEnabled(True)
                self.bmAxis.setCursor(Qt.CursorShape.ArrowCursor)
            else:
                return

    def reset_roi(self):
        if self.rect is not None:
            self.rect.remove()
            self.rect = None
            self.bmcanvas.draw()

        self.roixi = 0
        self.roizi = 0
        self.roixf = self.X*self.dx
        self.roizf = self.Z*self.dz
        self.programmatic_change = True
        self.roi_pbsetup()
        self.programmatic_change = False
        self.pbResetROI.setEnabled(False)
        self.bmax.set_xlim(0, self.X * self.dx)
        self.bmax.set_ylim(self.Z * self.dz, 0)
        self.update_frame()

    def roi_pbchanged(self):
        if self.programmatic_change:
            return
        else:
            if self.rect is not None:
                self.rect.remove()
                self.rect = None
                self.bmcanvas.draw()
            self.roixi = self.sbXi.value()
            self.roixf = self.sbXf.value()
            self.roizi = self.sbZi.value()
            self.roizf = self.sbZf.value()
            self.generate_roi()

    def roi_pbsetup(self):
        self.programmatic_change = True
        self.sbXi.setValue(0)
        self.sbXf.setValue(self.X * self.dx)
        self.sbXf.setMaximum(self.X * self.dx)
        self.sbZi.setValue(0)
        self.sbZf.setValue(self.Z * self.dz)
        self.sbZf.setMaximum(self.Z * self.dz)
        self.programmatic_change = False

    def roi_type(self, index):
        if self.rect is not None:
            self.rect.remove()
            self.rect = None
            self.bmcanvas.draw()
        self.generate_roi()

    def generate_roi(self):
        if self.cbROItype.currentText() == 'Plot Rectangle':

            self.rect = Rectangle((self.roixi, self.roizi),
                          self.roixf - self.roixi, self.roizf - self.roizi,
                          linewidth=2, edgecolor='red', facecolor='none', linestyle='dashed')
            self.bmax.set_xlim(0, self.X * self.dx)
            self.bmax.set_ylim(self.Z * self.dz,0)
            self.bmax.add_patch(self.rect)
            self.bmcanvas.draw()
        else:
            self.rect = Rectangle((self.roixi, self.roizi),
                          self.roixf - self.roixi, self.roizf - self.roizi,
                          linewidth=2, edgecolor='none', facecolor='none', linestyle='dashed')
            self.bmax.add_patch(self.rect)
            self.bmax.set_xlim(self.roixi, self.roixf)
            self.bmax.set_ylim(self.roizf, self.roizi)
            self.bmcanvas.draw()

# ----------------------------------------------BMode Processing--------------------------------------------------------

    def US_avg(self):
        if self.bmloaded == True:
            N = self.sbUSavg.value()
            aux = USavg(self.rfdata, N)
            self.rfdata = aux
            [self.Z, self.X, self.frames] = (self.rfdata.shape)
            self.sbFr.setValue(int(self.fr/2))
            self.rfdata_frames_update()
            self.update_frame()

    def US_iniavg(self):
        if self.bmloaded == True:
            N = self.sbUSiniavg.value()
            aux = USiniavg(self.rfdata, N)
            self.rfdata = aux
            [self.Z, self.X, self.frames] = (self.rfdata.shape)
            self.rfdata_frames_update()
            self.update_frame()

    def US_interp(self):
        if self.bmloaded == True:
            aux = USinterp(self.rfdata)
            self.rfdata = aux
            [self.Z, self.X, self.frames] = (self.rfdata.shape)
            self.dx = self.dx/2
            self.rfdata_frames_update()
            self.sbdx.setValue(self.dx)
            self.pbUSinterp.setEnabled(False)


# --------------------------------------------US Parameters Panel-------------------------------------------------------
    def sf_changed(self):
        self.fs = (self.sbSf.value() * 1000000)
        self.dz = (self.c / (2 * self.fs)) * 1000
        if self.bmloaded == True:
            self.reset_roi()

    def dx_changed(self):
        self.dx = self.sbdx.value()
        if self.bmloaded == True:
            self.reset_roi()

    def c_changed(self):
        self.c = self.sbC.value()
        self.dz = (self.c / (2 * self.fs)) * 1000
        if self.bmloaded == True:
            self.reset_roi()

    def fr_changed(self):
        self.fr = self.sbFr.value()
        self.dt = float(1/self.fr)
        self.sbLPfilter.setMaximum(int(self.fr/2))
        self.sbLPfilter.setValue(int(self.fr / 2))

    def fc_changed(self):
        self.fc = (self.sbFc.value() * 1000000)

# ---------------------------------------------Displacement Panel-------------------------------------------------------
    def set_dispsliders(self):
        # dimensions
        self.dispmax = np.amax(self.dispmap)
        self.dispmin = np.amin(self.dispmap)
        self.dispdz = (self.roizf - self.roizi) / self.dispZ
        self.dispdx = (self.roixf - self.roixi) / self.dispX

        self.slDisp.setMinimum(1)  # Set the minimum value
        self.slDisp.setMaximum(self.dispframes)
        #self.slDisp.setValue(1)

        self.slDispMax.setMinimum(0)  # Set the minimum value
        self.slDispMax.setMaximum(int(self.dispmax*100))
        self.slDispMax.setSingleStep(int(self.dispmax*10))
        self.slDispMax.setValue(int(self.dispmax*100))
        self.lbDispMax.setText(f'{self.dispmax: .2f}')

        self.slDispMin.setMinimum(int(self.dispmin*100))  # Set the minimum value
        self.slDispMin.setMaximum(0)
        self.slDispMin.setSingleStep(int(self.dispmin*10))
        self.slDispMin.setValue(int(self.dispmin*100))
        self.lbDispMin.setText(f'{self.dispmin: .2f}')

    def loupas(self):
        self.activate_guiBM(False)
        self.activate_guiDisp(False)
        particle = self.cbParticle.currentText()
        if particle == 'Velocity':
            velo = 1
            integrate = 0
            diff = 0
        elif particle == 'Displacement':
            velo = 0
            integrate = 0
            diff = 0
        elif particle == 'Integrate Velocity':
            velo = 1
            integrate = 1
            diff = 0
        elif particle == 'Diff Displacement':
            velo = 0
            integrate = 0
            diff = 1

        trkprm = {
            'fs': self.fs,
            'fc': self.fc,
            'c': self.c,
            'winsize': self.winsize,
            'overlap': self.overlap,
            'df': self.skipf,
            'ini': (self.inif - 1),
            'end': self.endf - 1,
            'velo': velo,
        }
        xipx = math.floor(self.roixi/self.dx)
        xfpx = math.floor(self.roixf / self.dx)
        zipx = math.floor(self.roizi / self.dz)
        zfpx = math.floor(self.roizf / self.dz)

        roi = {
            'xi': int(xipx),
            'xf': int(xfpx),
            'zi': int(zipx),
            'zf': int(zfpx),
        }

        IQbuffer = calculateIQ(self.rfdata, trkprm,roi)
        self.dispmap = framesLoupas(IQbuffer,trkprm)

        if integrate == 1:
            self.dispmap = integrateDisp(self.dispmap)

        if diff == 1:
            self.dispmap = diffDisp(self.dispmap)

        [self.dispZ,self.dispX,self.dispframes] = self.dispmap.shape

        self.set_dispsliders()
        self.activate_guiBM(True)
        self.activate_guiDisp(True)
        self.dispupdate_frame()

    def kasai(self):
        self.activate_guiBM(False)
        self.activate_guiDisp(False)
        particle = self.cbParticle.currentText()
        if particle == 'Velocity':
            velo = 1
        else:
            velo = 0
        trkprm = {
            'fs': self.fs,
            'fc': self.fc,
            'c': self.c,
            'winsize': self.winsize,
            'overlap': self.overlap,
            'df': self.skipf,
            'ini': (self.inif - 1),
            'end': self.endf-1,
            'velo': velo,
        }
        xipx = math.floor(self.roixi / self.dx)
        xfpx = math.floor(self.roixf / self.dx)
        zipx = math.floor(self.roizi / self.dz)
        zfpx = math.floor(self.roizf / self.dz)

        roi = {
            'xi': int(xipx),
            'xf': int(xfpx),
            'zi': int(zipx),
            'zf': int(zfpx),
        }

        IQbuffer = calculateIQ(self.rfdata, trkprm, roi)
        self.dispmap = framesKasai(IQbuffer, trkprm)

        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape

        self.set_dispsliders()
        self.activate_guiBM(True)
        self.activate_guiDisp(True)
        self.dispupdate_frame()

    def cross_cor(self):
        self.activate_guiBM(False)
        self.activate_guiDisp(False)
        particle = self.cbParticle.currentText()
        if particle == 'Velocity':
            velo = 1
        else:
            velo = 0
        trkprm = trkprm = {
            'fs': self.fs,
            'fc': self.fc,
            'c': self.c,
            'winsize': self.winsize,
            'overlap': self.overlap,
            'df': self.skipf,
            'ini': (self.inif - 1),
            'end': self.endf-1,
            'velo': velo,
        }
        xipx = math.floor(self.roixi / self.dx)
        xfpx = math.floor(self.roixf / self.dx)
        zipx = math.floor(self.roizi / self.dz)
        zfpx = math.floor(self.roizf / self.dz)

        roi = {
            'xi': int(xipx),
            'xf': int(xfpx),
            'zi': int(zipx),
            'zf': int(zfpx),
        }


        rfdata = cropRfdata(self.rfdata, trkprm, roi)
        self.dispmap = framesCrossCor(rfdata, trkprm)
        self.dispmap = self.dispmap*self.dz*1000
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape

        self.set_dispsliders()
        self.activate_guiBM(True)
        self.activate_guiDisp(True)
        self.dcDisp()
        self.dispupdate_frame()

    def tracking_prm(self):
        fdialog = Tracking_Dialog()

        fdialog.set_win(self.Z,self.winsize)
        fdialog.set_overlap(self.overlap)
        fdialog.set_skip(self.skipf)
        fdialog.set_end(self.frames)

        # Define variables to store ini and end values
        win = None
        overlap = None
        ini = None
        skip = None
        end = None

        # Connect the frames_index signal to update ini and end variables
        def update_tracking(w,o,i,s,e):
            nonlocal win,overlap,ini, skip, end
            win,overlap,ini, skip, end = (w,o,i,s,e)

        fdialog.trk_prm.connect(update_tracking)

        # Execute the dialog
        fdialog.exec()

        if (win == -1 and overlap == -1  and ini == -1  and skip == -1  and end == -1):
            return
        else:
            self.winsize = win
            self.overlap = overlap
            self.inif = ini
            self.skipf = skip
            self.endf = end

#------------------------------------------------Time Filtering---------------------------------------------------------
    def LPfilter(self):

        cutoff_freq = float(self.sbLPfilter.value());
        cutoff_freq = cutoff_freq/(float(self.fr/2))

        self.dispmap = DispLPfilter(self.dispmap,cutoff_freq)

        self.set_dispsliders()
        self.dispupdate_frame()

    def BPfilter(self):


        cutoff_freq1 = float(self.sbBPfilter1.value());
        cutoff_freq2 = float(self.sbLPfilter.value());
        if(cutoff_freq1 >= cutoff_freq2):
            cutoff_freq1 = cutoff_freq2-20
            self.sbBPfilter1.setValue(int(cutoff_freq1))
        cutoff_freq1 = cutoff_freq1 / (float(self.fr/2))
        cutoff_freq2 = cutoff_freq2 / (float(self.fr / 2))

        self.dispmap = DispBPfilter(self.dispmap,cutoff_freq1,cutoff_freq2)

        self.set_dispsliders()
        self.dispupdate_frame()

    def interpDisp(self):
        self.dispmap = Dispinterp(self.dispmap)
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.sbFr.setValue(int(self.fr * 2))
        self.set_dispsliders()
        self.dispupdate_frame()

    def frames_selection(self):
        self.dispmap = cropFrames(self.dispmap)
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.set_dispsliders()
        self.dispupdate_frame()

    def detrendDisp(self):
        self.dispmap = Dispdetrend(self.dispmap)
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.set_dispsliders()
        self.dispupdate_frame()

#-----------------------------------------------Spatial Filtering-------------------------------------------------------
    def medianDisp(self):
        self.dispmap = medianFilter(self.dispmap,self.mx,self.mz)
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.set_dispsliders()
        self.dispupdate_frame()

    def leeDisp(self):
        self.dispmap = leeFilter(self.dispmap,self.mx,self.mz)
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.set_dispsliders()
        self.dispupdate_frame()

    def frostDisp(self):

        hx = self.mx // 2
        hx_floor = np.floor(hx)
        hz = self.mx // 2
        hz_floor = np.floor(hz)

        # Create meshgrid-like arrays
        x_range = np.arange(-hz_floor, hz_floor + 1)
        z_range = np.arange(-hx_floor, hx_floor + 1)
        x, z = np.meshgrid(x_range, z_range)

        # Calculate matrix S
        S = np.sqrt(x ** 2 + z ** 2)


        # dll load
        dll_path = r".\C_functions\cmmus.dll"
        fdll = ctypes.CDLL(dll_path)
        fdll.frostfilter.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64),  # bfData
            np.ctypeslib.ndpointer(dtype=np.float64),  # bfData
            np.ctypeslib.ndpointer(dtype=np.float64),  # bfData
            np.ctypeslib.ndpointer(dtype=np.int32)  # prm
        ]

        prm = np.zeros((4))

        prm[0] = self.dispZ
        prm[1] = self.dispX
        prm[2] = self.mz
        prm[3] = self.mx
        dialog = QApplication.instance()
        progress_bar_window = ProgressBarWindow()
        progress_bar_window.center_on_screen()
        progress_bar_window.set_title("Frost Filter")
        progress_bar_window.show()
        inter = np.zeros((self.dispZ, self.dispX))
        interf = np.zeros((self.dispZ, self.dispX))
        for i in range(0, self.dispframes):
            inter = self.dispmap[:, :, i]
            fdll.frostfilter(interf, np.float64(inter.flatten()), np.float64(S.flatten()), np.int32(prm))
            progress_bar_window.set_progress_value(int((i / (self.frames)) * 100))
            progress_bar_window.repaint()
            dialog.processEvents()
            self.dispmap[:,:,i] = interf

        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.set_dispsliders()
        self.dispupdate_frame()

    def dircall(self):
        self.dispmap = directionalFilter(self.dispmap, self.dx, self.dz, [self.lc, self.hc], self.order, self.power, self.angle)
        [self.dispZ, self.dispX, self.dispframes] = self.dispmap.shape
        self.set_dispsliders()
        self.dispupdate_frame()

    def dcDisp(self):
        self.dispmap = removeDC(self.dispmap)
        self.set_dispsliders()
        self.dispupdate_frame()

    def filter_prm(self):
        fdialog = sfilter_Dialog()


        fdialog.set_mz(self.mz)
        fdialog.set_mx(self.mx)
        fdialog.set_angle(self.angle)
        fdialog.set_lc(self.lc)
        fdialog.set_hc(self.hc)
        fdialog.set_order(self.order)
        fdialog.set_power(self.power)

        # Define variables to store ini and end values
        mz = None
        mx = None
        angle = None
        lc = None
        hc = None
        order = None
        power = None

        # Connect the frames_index signal to update ini and end variables
        def update_filters(z,x,a,l,h,o,p):
            nonlocal mz,mx,angle,lc,hc,order,power
            mz,mx,angle,lc,hc,order,power = (z,x,a,l,h,o,p)

        fdialog.sfilter_prm.connect(update_filters)

        # Execute the dialog
        fdialog.exec()

        if (mz == -1 and mx == -1  and angle == -1  and lc == -1  and hc == -1  and order == -1  and power == -1):
            return
        else:
            self.mz = mz
            self.mx = mx
            self.angle = angle
            self.lc = lc
            self.hc = hc
            self.order = order
            self.power = power

# ------------------------------------------------Time Analysis---------------------------------------------------------
    def signal_point(self):
        if self.Signaltext is not None:
            self.Signaltext.remove()
            self.dispcanvas.draw()
        self.Signaltext = self.dispax.text(0.5, 0.5, "Select a Point",
                                    color='black', fontsize=12, ha='center', va='center', transform=self.dispax.transAxes, )
        self.dispcanvas.draw()
        self.dispAxis.setCursor(Qt.CursorShape.CrossCursor)
        self.connection_id2 = self.dispfig.canvas.mpl_connect('button_press_event', self.signal_selection)

    def signal_selection(self, event):
        if event.inaxes == self.dispax:
            if self.Signaltext is not None:
                self.Signaltext.remove()
                self.Signaltext = None
                self.dispcanvas.draw()
            if self.Signalmarker is not None:
                self.dispshow_image()
                self.graphax.clear()
                self.Signalmarker = None
            sgx, sgz = event.xdata, event.ydata

            sgx_px = int((sgx - self.roixi) / self.dispdx)
            sgz_px = int((sgz - self.roizi) / self.dispdz)
            self.Signalmarker = self.dispax.plot(sgx, sgz, 'ko')
            self.dispax.figure.canvas.draw()

            self.signal_time = np.arange(0, self.dispframes, 1)
            self.signal_time = self.signal_time * self.dt * 1000

            self.signal_amp = np.zeros((self.dispframes))
            self.signal_amp = np.squeeze(self.dispmap[sgz_px,sgx_px,:])

            self.graphfig.clf()
            self.graphax = self.graphfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])

            self.graphax.clear()
            self.graphfig.set_facecolor((0.94, 0.94, 0.94))
            self.graphax.set_facecolor((1, 1, 1))
            self.graphax.plot(self.signal_time, self.signal_amp, 'k-')
            self.graphax.set_ylabel('Amplitude [μm]', fontsize=9)
            self.graphax.set_xlabel('Time [ms]', fontsize=9)
            self.graphax.tick_params(axis='x', labelsize=9)
            self.graphax.tick_params(axis='y', labelsize=9)
            self.graphax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')

            self.graphcanvas.draw()

            self.dispfig.canvas.mpl_disconnect(self.connection_id2)
            self.dispAxis.setCursor(Qt.CursorShape.ArrowCursor)

            self.pbPSD.setEnabled(True)
            self.activate_velosliders(False)
            if self.velocursorflag == True:
                self.velo_cursor()

    def signal_PSD(self):

        [ps,freq] = powerSpectrum(self.signal_amp,self.fr)
        self.graphfig.clf()
        self.graphax = self.graphfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.graphax.clear()
        self.graphfig.set_facecolor((0.94, 0.94, 0.94))
        self.graphax.set_facecolor((1, 1, 1))
        self.graphax.plot(freq, ps, 'k-')
        self.graphax.set_ylabel('Amplitude [a.u.]', fontsize=9)
        self.graphax.set_xlabel('Frequency [Hz]', fontsize=9)
        self.graphax.tick_params(axis='x', labelsize=9)
        self.graphax.tick_params(axis='y', labelsize=9)
        self.graphax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')
        self.graphcanvas.draw()

        self.pbPSD.setEnabled(False)
        self.activate_velosliders(False)
        if self.velocursorflag == True:
            self.velo_cursor()

    def areaAnalysis(self):
        prm = {
            'dx': self.dispdx,
            'dz': self.dispdz,
            'dt': self.dt,
            'index': self.dispcurrent_frame,
            'roixi': self.roixi,
            'roixf': self.roixf,
            'roizi': self.roizi,
            'roizf': self.roizf,

        }
        areaDisp(self.dispmap, prm)

# ------------------------------------------------Group Velocity--------------------------------------------------------
    def ToF(self):
        cbmethod = self.cbGmethod.currentText()
        if cbmethod == 'Cross Correlation':
            gmethod = 'corr'
        elif cbmethod == 'Maximum':
            gmethod = 'max'
        elif cbmethod == 'Minimum':
            gmethod = 'min'
        dialog = QApplication.instance()
        progress_bar_window = ProgressBarWindow()
        progress_bar_window.center_on_screen()
        progress_bar_window.set_title("Group Velocity Calculation")
        progress_bar_window.show()
        #tlim = self.gtr
        N = self.N
        dist = self.dist
        tlim = self.dist / (self.gtr * 1000 * self.dt)
        [dZ,dX,dframes] = self.dispmap.shape
        self.roixi, self.roixf, self.roizf, self.roizi
        self.vx = np.zeros((dZ,dX))
        self.vz = np.zeros((dZ, dX))
        self.v = np.zeros((dZ, dX))
        dx = (self.roixf - self.roixi)/dX
        dz = (self.roizf - self.roizi)/dZ
        dimz = np.arange(0, dz * (dZ + 1), dz)
        dimx = np.arange(0, dx * (dX + 1), dx)
        stx = round(dist/dx)
        stz = round(dist/dz)
        wsizex = round((N - 1) / 2 * stx)
        wsizez = round((N - 1) / 2 * stz)
        for x in range(wsizex+1,dX-(wsizex+1)):
            for z in range(wsizez + 1, dZ - (wsizez + 1)):
                posx = np.zeros((N,N))
                posz = np.zeros((N, N))
                recorte = np.zeros((dframes,N,N))
                cx = 0
                for i in range(-wsizex,wsizex+1,stx):
                    cz = 0
                    for j in range(-wsizez, wsizez+1, stz):
                        recorte[:,cz,cx] = np.squeeze(self.dispmap[z + j,x + i,:])
                        posx[cz,cx] = dimx[x + i];
                        posz[cz,cx] = dimz[z + j];
                        cz = cz + 1;
                    cx = cx + 1;
                self.vx[z,x],self.vz[z,x],self.v[z,x] = swgvelo2(recorte, posx, posz, self.dt, tlim, gmethod)
            progress_bar_window.set_progress_value(int((x / (dX-(wsizex+1))) * 100))
            progress_bar_window.repaint()
            dialog.processEvents()

        progress_bar_window.close_window()
        self.velosetsliders()
        self.activate_velosliders(True)
        self.showVelo()
        self.velocursorflag = False

    def sToF(self,z,x):
        cbmethod = self.cbGmethod.currentText()
        if cbmethod == 'Cross Correlation':
            gmethod = 'corr'
        elif cbmethod == 'Maximum':
            gmethod = 'max'
        elif cbmethod == 'Minimum':
            gmethod = 'min'

        N = self.N
        dist = self.dist
        tlim = self.dist/(self.gtr*1000*self.dt)
        [dZ,dX,dframes] = self.dispmap.shape
        self.roixi, self.roixf, self.roizf, self.roizi
        dx = (self.roixf - self.roixi)/dX
        dz = (self.roizf - self.roizi)/dZ
        dimz = np.arange(self.roizi, self.roizf, dz)
        dimx = np.arange(self.roixi, self.roixf, dx)
        stx = round(dist/dx)
        stz = round(dist/dz)
        wsizex = round((N - 1) / 2 * stx)
        wsizez = round((N - 1) / 2 * stz)
        posx = np.zeros((N,N))
        posz = np.zeros((N, N))
        recorte = np.zeros((dframes,N,N))
        cx = 0
        for i in range(-wsizex,wsizex+1,stx):
            if(0 <= x + i <self.dispX):
                cz = 0
                for j in range(-wsizez, wsizez+1, stz):
                    if (0 <= z + j < self.dispZ):
                        recorte[:,cz,cx] = np.squeeze(self.dispmap[z + j,x + i,:])
                        posx[cz,cx] = dimx[x + i];
                        posz[cz,cx] = dimz[z + j];
                        self.dispax.plot(posx[cz,cx], posz[cz,cx], marker='x', linestyle='None', color='k', markersize=5)
                        cz = cz + 1;
                cx = cx + 1;
        self.dispax.figure.canvas.draw()
        t,sx,sz = swgvelo1(recorte, posx, posz, self.dt, tlim, gmethod)
        ht = np.argmin(np.abs(t - 0))
        s = np.zeros((len(sx)))
        for i in range(0,len(sx)):
            self.dispax.plot(sx[i], sz[i], marker='x', linestyle='None', color='r', markersize=5)
            s[i] = math.sqrt((sz[i] - sz[ht])*(sz[i] - sz[ht])  + (sx[i] - sx[ht])*(sx[i] - sx[ht]))
        self.dispax.plot(dimx[x], dimz[z], marker='x', linestyle='None', color='r', markersize=5)
        hs = np.argmin(np.abs(s - 0))  # Find index of element in s closest to 0
        # Multiply elements up to index hs (inclusive) by -1
        s[:hs + 1] *= -1
        self.sToF_graph(t,s)
        self.dispax.figure.canvas.draw()

    def sToF_point(self):
        self.dispshow_image()
        self.Signaltext = self.dispax.text(0.5, 0.5, "Select a Point",
                                    color='black', fontsize=12, ha='center', va='center', transform=self.dispax.transAxes, )
        self.dispcanvas.draw()
        self.dispAxis.setCursor(Qt.CursorShape.CrossCursor)
        self.connection_id2 = self.dispfig.canvas.mpl_connect('button_press_event', self.sToF_selection)

    def sToF_selection(self, event):
        if event.inaxes == self.dispax:
            if self.Signaltext is not None:
                self.Signaltext.remove()
                self.Signaltext = None
                #self.dispcanvas.draw()
            sgx, sgz = event.xdata, event.ydata

            sgx_px = int((sgx - self.roixi) / self.dispdx)
            sgz_px = int((sgz - self.roizi) / self.dispdz)

            #self.Signalmarker = self.dispax.plot(sgx, sgz, marker='x', linestyle='None', color='k', markersize=5)
            #self.dispax.figure.canvas.draw()

            self.dispfig.canvas.mpl_disconnect(self.connection_id2)
            self.dispAxis.setCursor(Qt.CursorShape.ArrowCursor)
            self.sToF(sgz_px,sgx_px)

    def sToF_graph(self,t,s):
        self.graphfig.clf()
        self.graphax = self.graphfig.add_axes(rect=[0.17, 0.14, 0.75, 0.75])
        self.graphax.clear()
        self.graphfig.set_facecolor((0.94, 0.94, 0.94))
        self.graphax.set_facecolor((1, 1, 1))
        self.graphax.plot(t*1000, s, marker='o', linestyle='None', color='k', markersize=5)
        self.graphax.set_ylabel('Distance [mm]', fontsize=9)
        self.graphax.set_xlabel('Time [ms]', fontsize=9)
        self.graphax.tick_params(axis='x', labelsize=9)
        self.graphax.tick_params(axis='y', labelsize=9)
        self.graphax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')

        if len(t) >= 3:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=np.RankWarning)
                coefficients = np.polyfit(t*1000, s, 1)
                v = coefficients[0]  # Coefficient of x^1 (slope)
                formatted_string = "{:.3f}".format(v)
                self.graphax.set_title('V = ' + formatted_string, fontsize=9)
                self.graphcanvas.draw()
                tfit = np.linspace(np.min(t), np.max(t), num=100)
                sfit = np.zeros((100))

                for i in range(0,100):
                    sfit[i] = 1000*tfit[i]*v
                self.graphax.plot(1000*tfit, sfit, 'r-')

        self.graphcanvas.draw()
        self.activate_velosliders(False)
        if self.velocursorflag == True:
            self.velo_cursor()

    def velosetsliders(self):
        self.programmatic_change = True
        self.slVeloMax.setMinimum(0)  # Set the minimum value
        self.slVeloMax.setMaximum(int(self.gtr * 100))
        self.slVeloMax.setSingleStep(int(10))
        self.slVeloMax.setValue(int(self.gtr * 100))
        self.lbVeloMax.setText(f'{self.gtr: .2f}')

        self.slVeloMin.setMinimum(-int(self.gtr * 100))  # Set the minimum value
        self.slVeloMin.setMaximum(0)
        self.slVeloMin.setSingleStep(int(10))
        self.slVeloMin.setValue(-int(self.gtr * 100))
        self.lbVeloMin.setText(f'{-self.gtr: .2f}')
        self.programmatic_change = False

    def showVelo(self):
        self.graphax.clear()
        self.graphfig.set_facecolor((0.94, 0.94, 0.94))
        self.graphax.set_facecolor((1, 1, 1))
        min = float(self.slVeloMin.value()) / 100
        max = float(self.slVeloMax.value()) / 100
        extent = [self.roixi, self.roixf, self.roizf, self.roizi]
        veloshow = self.cbVeloshow.currentText()
        if veloshow == 'Vx':
            self.graphax.imshow(self.vx, cmap='jet',
                           extent=extent, vmin=min, vmax=max, interpolation='none')
        elif veloshow == 'Vz':
            self.graphax.imshow(self.vz, cmap='jet',
                           extent=extent, vmin=min, vmax=max, interpolation='none')
        elif veloshow == 'Velocity':
            self.graphax.imshow(self.v, cmap='jet',
                           extent=extent, vmin=min, vmax=max, interpolation='none')

        self.graphax.set_xlabel('Lateral [mm]', fontsize=9)
        self.graphax.set_ylabel('Axial [mm]', fontsize=9)
        self.graphax.tick_params(axis='x', labelsize=9)
        self.graphax.tick_params(axis='y', labelsize=9)
        self.graphax.set_aspect('equal',adjustable = 'box')
        self.graphcanvas.draw()
        #self.graphax.set_aspect('equal')

    def veloupdate_frame(self):
        min = float(self.slVeloMin.value()) / 100
        max = float(self.slVeloMax.value()) / 100
        self.lbVeloMax.setText(f'{max: .2f}')
        self.lbVeloMin.setText(f'{min: .2f}')
        self.showVelo()

    def veloslider_changed(self):
        if self.programmatic_change == False:
            self.veloupdate_frame()
        else:
            return

    def gvelo_set(self):
        fdialog = Gvelocity_Dialog()

        fdialog.set_dist(self.dist)
        fdialog.set_kernel(self.N)
        fdialog.set_tr(self.gtr)

        # Define variables to store ini and end values
        dist = None
        N = None
        tr = None

        # Connect the frames_index signal to update ini and end variables
        def update_tracking(d,n,t):
            nonlocal dist, N, tr
            dist, N, tr = (d,n,t)

        fdialog.gvelo_prm.connect(update_tracking)

        # Execute the dialog
        fdialog.exec()

        if (dist == -1 and N == -1  and tr == -1):
            return
        else:
            self.dist = dist
            self.gtr = tr
            if N % 2 == 0:
                self.N = int(N+1)
            else:
                self.N = int(N)

        self.velosetsliders()

    def veloshow_change(self):
        self.showVelo()

    def on_motion_velo_cursor(self,event):
        if event.xdata is not None and event.ydata is not None:
            x = event.xdata
            y = event.ydata

            veloshow = self.cbVeloshow.currentText()
            if veloshow == 'Vx':
                x_data = int((x - self.roixi) / (self.roixf - self.roixi) * self.vx.shape[1])
                y_data = int(self.vx.shape[0] - (y - self.roizf) / (self.roizi - self.roizf) * self.vx.shape[0])
                if 0 <= x < self.vx.shape[1] and 0 <= y < self.vx.shape[0]:
                    value = self.vx[y_data, x_data]
                    self.lbVeloval.setText(f'V = {value: .2f} m/s')
                    self.lbVeloval.setAlignment(Qt.AlignmentFlag.AlignHCenter)
            elif veloshow == 'Vz':
                x_data = int((x - self.roixi) / (self.roixf - self.roixi) * self.vz.shape[1])
                y_data = int(self.vz.shape[0] - (y - self.roizf) / (self.roizi - self.roizf) * self.vz.shape[0])
                if 0 <= x < self.vz.shape[1] and 0 <= y < self.vz.shape[0]:
                    value = self.vz[y_data, x_data]
                    self.lbVeloval.setText(f'V = {value: .2f} m/s')
                    self.lbVeloval.setAlignment(Qt.AlignmentFlag.AlignHCenter)
            elif veloshow == 'Velocity':
                x_data = int((x - self.roixi) / (self.roixf - self.roixi) * self.v.shape[1])
                y_data = int(self.v.shape[0] - (y - self.roizf) / (self.roizi - self.roizf) * self.v.shape[0])
                if 0 <= x < self.v.shape[1] and 0 <= y < self.v.shape[0]:
                    value = self.v[y_data, x_data]
                    self.lbVeloval.setText(f'V = {value: .2f} m/s')
                    self.lbVeloval.setAlignment(Qt.AlignmentFlag.AlignHCenter)

    def velo_cursor(self):
        if self.velocursorflag:
            self.graphfig.canvas.mpl_disconnect(self.connection_velo)
            self.graphAxis.setCursor(Qt.CursorShape.ArrowCursor)
            self.lbVeloval.setVisible(False)
            self.pbVelocursor.setText('Cursor off')
            self.velocursorflag = False
        else:
            self.connection_velo = self.graphfig.canvas.mpl_connect('motion_notify_event', self.on_motion_velo_cursor)
            self.graphAxis.setCursor(Qt.CursorShape.CrossCursor)
            self.lbVeloval.setVisible(True)
            self.pbVelocursor.setText('Cursor on')
            self.velocursorflag = True

    def veloROI(self):
        prm = {
            'dx': self.dispdx,
            'dz': self.dispdz,
            'dt': self.dt,
        }
        veloshow = self.cbVeloshow.currentText()
        if veloshow == 'Vx':
            areaVelo(self.vx, prm)
        elif veloshow == 'Vz':
            areaVelo(self.vz, prm)
        elif veloshow == 'Velocity':
            areaVelo(self.v, prm)

# ------------------------------------------------Group Velocity--------------------------------------------------------

    def phaseVelo1(self):
        prm = {
            'dx': self.dispdx,
            'dz': self.dispdz,
            'dt': self.dt,
            'index': self.dispcurrent_frame,
            'zeropad': self.zeropad,
            'tr': self.pvelotr,
            'freqcut': self.pvelofc,
            'df': self.pvelodf,
            'density': self.pvelodensity,
         }
        [auxfreqdata, auxdispersion] = phaseVelo(self.dispmap, prm)
        if (len(auxfreqdata)==1 and len(auxdispersion)==1):
            return
        else:
            # Initialize an array to store the closest values
            target_multiples = np.arange(0, np.max(auxfreqdata), self.pvelodf)
            self.freqdata = []
            self.dispersion = []
        

            # Iterate over each target multiple
            for multiple in target_multiples:
                # Find the index of the closest value in the original array to the target multiple
                closest_index = np.argmin(np.abs(auxfreqdata - multiple))
                # Add the closest value to the list
                self.freqdata.append(auxfreqdata[closest_index])
                self.dispersion.append(auxdispersion[closest_index])
            self.dispersion_show()

    def dispersion_show(self):

        self.graphfig.clf()
        self.graphax = self.graphfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.graphax.clear()
        self.graphfig.set_facecolor((0.94, 0.94, 0.94))
        self.graphax.set_facecolor((1, 1, 1))
        self.graphax.plot(self.freqdata,self.dispersion, 'ko', markersize=2)
        self.graphax.set_ylabel('Velocity [m/s]', fontsize=9)
        self.graphax.set_xlabel('Frequency [Hz]', fontsize=9)
        self.graphax.tick_params(axis='x', labelsize=9)
        self.graphax.tick_params(axis='y', labelsize=9)
        self.graphax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')
        self.graphcanvas.draw()

    def pvelo_prm(self):
        fdialog = Pvelocity_Dialog()

        fdialog.set_zeropad(self.zeropad)
        fdialog.set_tr(self.pvelotr)
        fdialog.set_freqcut(self.pvelofc, int((1/self.dt)/2))
        fdialog.set_df(self.pvelodf)
        fdialog.set_density(self.pvelodensity)

        # Define variables to store ini and end values
        zeropad = None
        tr = None
        freqcut = None
        df = None
        density = None

        # Connect the frames_index signal to update ini and end variables
        def update_pvelo(z, t, f, d, den):
            nonlocal zeropad, tr, freqcut, df, density
            zeropad, tr, freqcut, df, density = (z, t, f, d, den)

        fdialog.pvelo_prm.connect(update_pvelo)

        # Execute the dialog
        fdialog.exec()

        if (zeropad == -1 and tr == -1  and freqcut == -1  and df == -1  and density == -1):
            return
        else:
            self.zeropad = zeropad
            self.pvelotr = tr
            self.pvelofc = freqcut
            self.pvelodf = df
            self.pvelodensity = density

app = QtWidgets.QApplication(sys.argv)
window = pyMMUS()
window.show()
app.exec()