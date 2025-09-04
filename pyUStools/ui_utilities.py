import numpy
import numpy as np
from PyQt6.QtWidgets import *
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtGui import QGuiApplication, QFont
from PyQt6.QtCore import pyqtSignal, Qt

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import Rectangle, RectangleSelector

class ProgressBarWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        vbox = QVBoxLayout()

        # Create a QProgressBar
        self.progress_bar = QProgressBar(self)
        vbox.addWidget(self.progress_bar)

        self.setLayout(vbox)

        self.setGeometry(300, 500, 300, 50)

        self.setWindowFlags(self.windowFlags() & ~Qt.WindowType.WindowCloseButtonHint & ~Qt.WindowType.WindowMaximizeButtonHint & ~Qt.WindowType.WindowMinimizeButtonHint)
    def center_on_screen(self):
        # Get the geometry of the screen
        primary_screen = QGuiApplication.primaryScreen()
        screen_rect = primary_screen.availableGeometry()

        # Center the window on the screen
        self.move((screen_rect.width() - self.width()) // 2, (screen_rect.height() - self.height()) // 2)

    def set_title(self, title):
        self.setWindowTitle(title)
    def set_progress_value(self, value):
        self.progress_bar.setValue(value)

    def close_window(self):
        self.close()

class Frames_Dialog(QDialog):
    frames_index = pyqtSignal(int, int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Frames Select")
        self.resize(300, 200)  # Adjust size
        self.center_on_screen()

        self.lbIni = QLabel("Initial Frame:", self)
        self.lbIni.setGeometry(20, 20, 100, 30)  # Adjust position and size

        self.sbIni = QSpinBox(self)
        self.sbIni.setGeometry(130, 20, 100, 30)  # Adjust position and size
        self.sbIni.setMinimum(0)  # Set minimum value
        self.sbIni.setMaximum(100)  # Set maximum value

        self.lbEnd = QLabel("Final Frame:", self)
        self.lbEnd.setGeometry(20, 70, 100, 30)  # Adjust position and size

        self.sbEnd = QSpinBox(self)
        self.sbEnd.setGeometry(130, 70, 100, 30)  # Adjust position and size
        self.sbEnd.setMinimum(0)  # Set minimum value
        self.sbEnd.setMaximum(100)  # Set maximum value

        self.pbOK = QPushButton("OK", self)
        self.pbOK.setGeometry(50, 120, 75, 30)  # Adjust position and size
        self.pbOK.clicked.connect(self.get_values)

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(175, 120, 75, 30)  # Adjust position and size
        self.pbCancel.clicked.connect(self.cancel_action)


    def center_on_screen(self):
        # Get the geometry of the screen
        primary_screen = QGuiApplication.primaryScreen()
        screen_rect = primary_screen.availableGeometry()

        # Center the window on the screen
        self.move((screen_rect.width() - self.width()) // 2, (screen_rect.height() - self.height()) // 2)

    def get_values(self):
        ini = self.sbIni.value()
        end = self.sbEnd.value()
        self.frames_index.emit(ini, end)
        self.accept()
    def cancel_action(self):
        self.frames_index.emit(-1, -1)
        self.accept()

class SaveDisp_Dialog(QDialog):
    saveas = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Save disp as:")
        self.resize(300, 120)  # Adjust size

        self.pbFigure = QPushButton("Figure", self)
        self.pbFigure.setGeometry(20, 30, 75, 30)
        self.pbFigure.clicked.connect(lambda: self.choice_clicked("Figure"))

        self.pbMovie = QPushButton("Movie", self)
        self.pbMovie.setGeometry(100, 30, 75, 30)
        self.pbMovie.clicked.connect(lambda: self.choice_clicked("Movie"))

        self.pbData = QPushButton("Data", self)
        self.pbData.setGeometry(180, 30, 75, 30)
        self.pbData.clicked.connect(lambda: self.choice_clicked("Data"))

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(180, 80, 75, 30)
        self.pbCancel.clicked.connect(lambda: self.choice_clicked("Cancel"))

    def choice_clicked(self, choice):
        self.saveas.emit(choice)
        self.accept()  # Close the dialog

class SaveVelo_Dialog(QDialog):
    saveas = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Save disp as:")
        self.resize(300, 120)  # Adjust size

        self.pbFigure = QPushButton("Figure", self)
        self.pbFigure.setGeometry(20, 30, 75, 30)
        self.pbFigure.clicked.connect(lambda: self.choice_clicked("Figure"))

        self.pbData = QPushButton("Data", self)
        self.pbData.setGeometry(100, 30, 75, 30)
        self.pbData.clicked.connect(lambda: self.choice_clicked("Data"))

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(180, 80, 75, 30)
        self.pbCancel.clicked.connect(lambda: self.choice_clicked("Cancel"))

    def choice_clicked(self, choice):
        self.saveas.emit(choice)
        self.accept()  # Close the dialog

class Tracking_Dialog(QDialog):
    trk_prm = pyqtSignal(int, int, int, int, int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Tracking Parameters")
        self.resize(300, 200)  # Adjust size
        self.center_on_screen()

        #Tracking
        self.lbwinsize = QLabel("Win Size [px]:", self)
        self.lbwinsize.setGeometry(0, 12, 81, 16)  # Adjust position and size
        self.lbwinsize.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbwin = QSpinBox(self)
        self.sbwin.setGeometry(90, 10, 61, 22)  # Adjust position and size
        self.sbwin.setMinimum(5)  # Set minimum value


        self.lboverlap = QLabel("Overlap [%]:", self)
        self.lboverlap.setGeometry(0, 52, 81, 16)  # Adjust position and size
        self.lboverlap.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sboverlap = QSpinBox(self)
        self.sboverlap.setGeometry(90, 50, 61, 22)  # Adjust position and size
        self.sboverlap.setMinimum(5)  # Set minimum value
        self.sboverlap.setMaximum(95)  # Set maximum value


        #frames
        self.lbIni = QLabel("Ini:", self)
        self.lbIni.setGeometry(140, 12, 61, 16)  # Adjust position and size
        self.lbIni.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbIni = QSpinBox(self)
        self.sbIni.setGeometry(210, 10, 61, 22)  # Adjust position and size
        self.sbIni.setMinimum(1)  # Set minimum value
        self.sbIni.setValue(2)

        self.lbSkip = QLabel("Step:", self)
        self.lbSkip.setGeometry(140, 52, 61, 16)  # Adjust position and size
        self.lbSkip.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbSkip = QSpinBox(self)
        self.sbSkip.setGeometry(210, 50, 61, 22)  # Adjust position and size
        self.sbSkip.setMinimum(1)  # Set minimum value
        self.sbSkip.setMaximum(10)  # Set maximum value


        self.lbEnd = QLabel("End:", self)
        self.lbEnd.setGeometry(140, 92, 61, 16)  # Adjust position and size
        self.lbEnd.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbEnd = QSpinBox(self)
        self.sbEnd.setGeometry(210, 90, 61, 22)  # Adjust position and size
        self.sbEnd.setMinimum(1)  # Set minimum value

        #Ok, Cancel
        self.pbOK = QPushButton("OK", self)
        self.pbOK.setGeometry(110, 140, 75, 30)  # Adjust position and size
        self.pbOK.clicked.connect(self.get_values)

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(200, 140, 75, 30)  # Adjust position and size
        self.pbCancel.clicked.connect(self.cancel_action)

    def center_on_screen(self):
        # Get the geometry of the screen
        primary_screen = QGuiApplication.primaryScreen()
        screen_rect = primary_screen.availableGeometry()

        # Center the window on the screen
        self.move((screen_rect.width() - self.width()) // 2, (screen_rect.height() - self.height()) // 2)

    def get_values(self):
        winsize = self.sbwin.value()
        overlap = self.sboverlap.value()
        ini = self.sbIni.value()
        skip = self.sbSkip.value()
        end = self.sbEnd.value()
        self.trk_prm.emit(winsize, overlap, ini, skip, end)
        self.accept()

    def cancel_action(self):
        self.trk_prm.emit(-1, -1,-1,-1,-1)
        self.accept()

    def set_win(self,max,value):
        self.sbwin.setValue(value)
        self.sbwin.setMaximum(max)

    def set_overlap(self,value):
        self.sboverlap.setValue(value)

    def set_skip(self,value):
        self.sbSkip.setValue(value)

    def set_end(self,value):
        self.sbEnd.setMaximum(value)  # Set maximum value
        self.sbEnd.setValue(value)
        self.sbIni.setMaximum(value-1)  # Set maximum value

class sfilter_Dialog(QDialog):
    sfilter_prm = pyqtSignal(int, int, int, float, float, int, float)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Spatial Filter Parameters")
        self.resize(320, 230)  # Adjust size
        self.center_on_screen()

        #General Filters

        self.lbGenfilt = QLabel("General Filters:", self)
        self.lbGenfilt.setGeometry(40, 8, 91, 16)
        self.lbGenfilt.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.lbMz = QLabel("Mz [px]:", self)
        self.lbMz.setGeometry(0, 30, 61, 16)
        self.lbMz.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbMz = QSpinBox(self)
        self.sbMz.setGeometry(70, 28, 61, 22)
        self.sbMz.setMinimum(1)
        self.sbMz.setSingleStep(2)
        self.sbMz.setMaximum(31)

        self.lbMx = QLabel("Mx [px]:", self)
        self.lbMx.setGeometry(140, 30, 61, 16)
        self.lbMx.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbMx = QSpinBox(self)
        self.sbMx.setGeometry(210, 28, 61, 22)
        self.sbMx.setMinimum(1)
        self.sbMx.setSingleStep(2)
        self.sbMx.setMaximum(31)



        #Directional Filter

        self.lbDirfilt = QLabel("Directional Filter:", self)
        self.lbDirfilt.setGeometry(40, 68, 91, 16)
        self.lbDirfilt.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.lbAngle = QLabel("Angle [º]:", self)
        self.lbAngle.setGeometry(0, 92, 61, 16)
        self.lbAngle.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbAngle = QSpinBox(self)
        self.sbAngle.setGeometry(70, 90, 61, 22)
        self.sbAngle.setMinimum(0)
        self.sbAngle.setSingleStep(1)
        self.sbAngle.setMaximum(360)


        self.lbLc = QLabel("Low cut:", self)
        self.lbLc.setGeometry(0, 122, 61, 16)
        self.lbLc.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbLc = QDoubleSpinBox(self)
        self.sbLc.setGeometry(70, 120, 61, 22)
        self.sbLc.setMinimum(0.01)
        self.sbLc.setSingleStep(0.01)
        self.sbLc.setMaximum(3)


        self.lbHc = QLabel("High cut:", self)
        self.lbHc.setGeometry(140, 122, 61, 16)
        self.lbHc.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbHc = QDoubleSpinBox(self)
        self.sbHc.setGeometry(210, 120, 61, 22)
        self.sbHc.setMinimum(0.01)
        self.sbHc.setSingleStep(0.1)
        self.sbHc.setMaximum(4)


        self.lbOrder = QLabel("Order:", self)
        self.lbOrder.setGeometry(0, 152, 61, 16)
        self.lbOrder.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbOrder = QSpinBox(self)
        self.sbOrder.setGeometry(70, 150, 61, 22)
        self.sbOrder.setMinimum(1)
        self.sbOrder.setSingleStep(1)
        self.sbOrder.setMaximum(10)


        self.lbPower = QLabel("Power:", self)
        self.lbPower.setGeometry(140, 152, 61, 16)
        self.lbPower.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbPower = QDoubleSpinBox(self)
        self.sbPower.setGeometry(210, 150, 61, 22)
        self.sbPower.setMinimum(0.01)
        self.sbPower.setSingleStep(0.01)
        self.sbPower.setMaximum(2)


        #Ok, Cancel
        self.pbOK = QPushButton("OK", self)
        self.pbOK.setGeometry(110, 180, 75, 30)
        self.pbOK.clicked.connect(self.get_values)

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(200, 180, 75, 30)
        self.pbCancel.clicked.connect(self.cancel_action)


    def center_on_screen(self):
        # Get the geometry of the screen
        primary_screen = QGuiApplication.primaryScreen()
        screen_rect = primary_screen.availableGeometry()

        # Center the window on the screen
        self.move((screen_rect.width() - self.width()) // 2, (screen_rect.height() - self.height()) // 2)

    def get_values(self):
        mz = self.sbMz.value()
        mx = self.sbMx.value()
        angle = self.sbAngle.value()
        lc = self.sbLc.value()
        hc = self.sbHc.value()
        order = self.sbOrder.value()
        power = self.sbPower.value()
        self.sfilter_prm.emit(mz, mx, angle, lc, hc, order, power)
        self.accept()

    def cancel_action(self):
        self.sfilter_prm.emit(-1, -1,-1,-1,-1, -1,-1)
        self.accept()

    def set_mz(self,value):
        self.sbMz.setValue(value)

    def set_mx(self,value):
        self.sbMx.setValue(value)

    def set_lc(self,value):
        self.sbLc.setValue(value)

    def set_hc(self,value):
        self.sbHc.setValue(value)

    def set_angle(self,value):
        self.sbAngle.setValue(value)

    def set_order(self,value):
        self.sbOrder.setValue(value)

    def set_power(self,value):
        self.sbPower.setValue(value)

class Gvelocity_Dialog(QDialog):
    gvelo_prm = pyqtSignal(float, int, float)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Group Velocity Parameters")
        self.resize(300, 200)  # Adjust size
        self.center_on_screen()

        #Tracking
        self.lbdist = QLabel("Distance [mm]:", self)
        self.lbdist.setGeometry(40, 12, 91, 16)  # Adjust position and size
        self.lbdist.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbdist = QDoubleSpinBox(self)
        self.sbdist.setGeometry(150, 10, 61, 22)  # Adjust position and size
        self.sbdist.setMinimum(0.1)  # Set minimum value
        self.sbdist.setMaximum(10)
        self.sbdist.setSingleStep(0.1)


        self.lbkernel = QLabel("Kernel [px]:", self)
        self.lbkernel.setGeometry(70, 52, 61, 16)  # Adjust position and size
        self.lbkernel.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbkernel = QSpinBox(self)
        self.sbkernel.setGeometry(150, 50, 61, 22)  # Adjust position and size
        self.sbkernel.setMinimum(3)  # Set minimum value
        self.sbkernel.setMaximum(31)  # Set maximum value
        self.sbkernel.setSingleStep(2)


        #frames
        self.lbTr = QLabel("Threshold [m/s]:", self)
        self.lbTr.setGeometry(40, 92, 91, 16)  # Adjust position and size
        self.lbTr.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbTr = QDoubleSpinBox(self)
        self.sbTr.setGeometry(150, 90, 61, 22)  # Adjust position and size
        self.sbTr.setMinimum(0)  # Set minimum value
        self.sbTr.setMaximum(30)
        self.sbTr.setSingleStep(0.5)

        #Ok, Cancel
        self.pbOK = QPushButton("OK", self)
        self.pbOK.setGeometry(110, 140, 75, 30)  # Adjust position and size
        self.pbOK.clicked.connect(self.get_values)

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(200, 140, 75, 30)  # Adjust position and size
        self.pbCancel.clicked.connect(self.cancel_action)


    def center_on_screen(self):
        # Get the geometry of the screen
        primary_screen = QGuiApplication.primaryScreen()
        screen_rect = primary_screen.availableGeometry()

        # Center the window on the screen
        self.move((screen_rect.width() - self.width()) // 2, (screen_rect.height() - self.height()) // 2)

    def get_values(self):
        dist = self.sbdist.value()
        N = self.sbkernel.value()
        tr = self.sbTr.value()
        self.gvelo_prm.emit(dist, N, tr)
        self.accept()
    def cancel_action(self):
        self.gvelo_prm.emit(-1, -1,-1)
        self.accept()

    def set_dist(self,value):
        self.sbdist.setValue(value)


    def set_kernel(self,value):
        self.sbkernel.setValue(value)

    def set_tr(self,value):
        self.sbTr.setValue(value)

class Pvelocity_Dialog(QDialog):
    pvelo_prm = pyqtSignal(int, int, int, int, int)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Phase Velocity Parameters")
        self.resize(400, 200)  # Adjust size
        self.center_on_screen()

        #zeropad
        self.lbzeropad = QLabel("Zero Pad [px]:", self)
        self.lbzeropad.setGeometry(20, 22, 91, 20)  # Adjust position and size
        self.lbzeropad.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbzeropad = QSpinBox(self)
        self.sbzeropad.setGeometry(120, 20, 61, 22)  # Adjust position and size
        self.sbzeropad.setMinimum(5)
        self.sbzeropad.setMaximum(5000)

        #threshold
        self.lbtr = QLabel("Threshold [dB]:", self)
        self.lbtr.setGeometry(10, 62, 101, 20)  # Adjust position and size
        self.lbtr.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbtr = QSpinBox(self)
        self.sbtr.setGeometry(120, 60, 61, 22)  # Adjust position and size
        self.sbtr.setMinimum(-100)  # Set minimum value
        self.sbtr.setMaximum(0)  # Set maximum value

        #freq cut
        self.lbfreqcut = QLabel("Freq Cut [Hz]:", self)
        self.lbfreqcut.setGeometry(0, 102, 111, 20)  # Adjust position and size
        self.lbfreqcut.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbfreqcut = QSpinBox(self)
        self.sbfreqcut.setGeometry(120, 100, 61, 22)  # Adjust position and size
        self.sbfreqcut.setMinimum(10)  # Set minimum value

        #df
        self.lbdf = QLabel("df [Hz]:", self)
        self.lbdf.setGeometry(230, 22, 61, 20)  # Adjust position and size
        self.lbdf.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbdf = QSpinBox(self)
        self.sbdf.setGeometry(300, 20, 61, 22)  # Adjust position and size
        self.sbdf.setMinimum(1)  # Set minimum value
        self.sbdf.setMaximum(100)  # Set maximum value

        #density
        self.lbdensity = QLabel("Density [Kg/m3]:", self)
        self.lbdensity.setGeometry(200, 62, 91, 20)  # Adjust position and size
        self.lbdensity.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.sbdensity = QSpinBox(self)
        self.sbdensity.setGeometry(300, 60, 61, 22)  # Adjust position and size
        self.sbdensity.setMinimum(100)  # Set minimum value
        self.sbdensity.setMaximum(10000)


        #Ok, Cancel
        self.pbOK = QPushButton("OK", self)
        self.pbOK.setGeometry(200, 140, 75, 23)  # Adjust position and size
        self.pbOK.clicked.connect(self.get_values)

        self.pbCancel = QPushButton("Cancel", self)
        self.pbCancel.setGeometry(286, 140, 75, 23)  # Adjust position and size
        self.pbCancel.clicked.connect(self.cancel_action)

    def center_on_screen(self):
        # Get the geometry of the screen
        primary_screen = QGuiApplication.primaryScreen()
        screen_rect = primary_screen.availableGeometry()

        # Center the window on the screen
        self.move((screen_rect.width() - self.width()) // 2, (screen_rect.height() - self.height()) // 2)

    def get_values(self):
        zeropad = self.sbzeropad.value()
        tr = self.sbtr.value()
        freqcut = self.sbfreqcut.value()
        df = self.sbdf.value()
        density = self.sbdensity.value()
        self.pvelo_prm.emit(zeropad, tr, freqcut, df, density)
        self.accept()

    def cancel_action(self):
        self.pvelo_prm.emit(-1, -1,-1,-1,-1)
        self.accept()

    def set_zeropad(self,value):
        self.sbzeropad.setValue(value)

    def set_tr(self,value):
        self.sbtr.setValue(value)

    def set_freqcut(self,value,mvalue):
        self.sbfreqcut.setMaximum(mvalue)
        self.sbfreqcut.setValue(value)

    def set_df(self,value):
        self.sbdf.setValue(value)

    def set_density(self,value):
        self.sbdensity.setValue(value)



class Disp_area(QDialog):

    def __init__(self, dispmap, prm, parent=None):
        super().__init__()
        self.initUI(dispmap,prm)

    def initUI(self, dispmap, prm):
        self.dispmap = dispmap
        [self.Z, self.X, self.frames] = (self.dispmap.shape)


        self.dx = prm['dx']
        self.dz = prm['dz']
        self.dt = prm['dt']
        self.index = prm['index']

        self.roixi = prm['roixi']
        self.roizi = prm['roizi']
        self.roixf = prm['roixf']
        self.roizf = prm['roizf']

        self.imroixi = prm['roixi']
        self.imroizi = prm['roizi']
        self.imroixf = prm['roixf']
        self.imroizf = prm['roizf']

        self.setWindowTitle("Area Analyzer")
        self.resize(1000, 700)

        self.max = np.amax(self.dispmap)
        self.min = np.amin(self.dispmap)

        #Disp Axis
        self.dispAxis = QWidget(self)
        self.dispAxis.setObjectName("dispAxis")
        self.dispAxis.setGeometry(30, 10, 591, 611)
        layout = QVBoxLayout(self.dispAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.dispfig = plt.Figure()
        self.dispfig.set_facecolor((0.94, 0.94, 0.94))
        self.dispax = self.dispfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.dispax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispax.tick_params(axis='x', labelsize=9)
        self.dispax.tick_params(axis='y', labelsize=9)
        self.dispcanvas = FigureCanvas(self.dispfig)
        layout.addWidget(self.dispcanvas)
        layout.addStretch()

        #Disp Slider
        self.slDisp = QSlider(self)
        self.slDisp.setGeometry(50, 620, 501, 16)
        self.slDisp.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.slDisp.setObjectName("slDisp")
        self.slDisp.setEnabled(False)
        self.programmatic_change = True
        self.slDisp.setMinimum(1)  # Set the minimum value
        self.slDisp.setMaximum(self.frames)
        self.slDisp.setValue(self.index)
        self.programmatic_change = False
        self.slDisp.valueChanged.connect(self.slider_changed)

        self.lbDisp = QLabel(self)
        self.lbDisp.setGeometry(560, 620, 61, 20)
        self.lbDisp.setObjectName("lbDisp")
        self.lbDisp.setText(f'{self.index} / {self.frames} ')

        self.pbResetROI = QPushButton(self)
        self.pbResetROI.setGeometry(50, 640, 75, 23)
        self.pbResetROI.setObjectName("pbResetROI")
        self.pbResetROI.setText("Reset ROI")
        self.pbResetROI.clicked.connect(self.select_roi)

        # Signal Axis
        self.signalAxis = QWidget(self)
        self.signalAxis.setObjectName("Widget3")
        self.signalAxis.setGeometry(640, 10, 331, 311)
        layout = QVBoxLayout(self.signalAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.signalfig = plt.Figure()
        self.signalfig.set_facecolor((0.94, 0.94, 0.94))
        self.signalax = self.signalfig.add_axes(rect=[0.22, 0.14, 0.75, 0.7])
        self.signalax.set_facecolor((0, 0, 0))
        self.signalax.set_xlabel('Time [ms]', fontsize=9)
        self.signalax.set_ylabel('Displacement [μm]', fontsize=9)
        self.signalax.tick_params(axis='x', labelsize=9)
        self.signalax.tick_params(axis='y', labelsize=9)
        self.signalcanvas = FigureCanvas(self.signalfig)
        layout.addWidget(self.signalcanvas)
        layout.addStretch()


        #PSD Axis
        self.psAxis = QWidget(self)
        self.psAxis.setObjectName("psAxis")
        self.psAxis.setGeometry(640, 330, 331, 311)
        layout = QVBoxLayout(self.psAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.psfig = plt.Figure()
        self.psfig.set_facecolor((0.94, 0.94, 0.94))
        self.psax = self.psfig.add_axes(rect=[0.22, 0.14, 0.75, 0.7])
        self.psax.set_facecolor((0, 0, 0))
        self.psax.set_xlabel('Frequency [Hz]', fontsize=9)
        self.psax.set_ylabel('Amplitude [a.u.]', fontsize=9)
        self.psax.tick_params(axis='x', labelsize=9)
        self.psax.tick_params(axis='y', labelsize=9)
        self.pscanvas = FigureCanvas(self.psfig)
        layout.addWidget(self.pscanvas)
        layout.addStretch()

        self.slDispMax = QSlider(self)
        self.slDispMax.setObjectName("slDispMax")
        self.slDispMax.setGeometry(577, 160, 22, 105)
        self.slDispMin = QSlider(self)
        self.slDispMin.setObjectName("slDispMin")
        self.slDispMin.setGeometry(577, 270, 22, 105)

        self.lbDispMax = QLabel(self)
        self.lbDispMax.setGeometry(600, 150, 61, 20)
        self.lbDispMax.setObjectName("lbDispMax")
        self.lbDispMin = QLabel(self)
        self.lbDispMin.setGeometry(600, 364, 61, 20)
        self.lbDispMin.setObjectName("lbDispMin")

        self.slDispMax.valueChanged.connect(self.dispslider_changed)
        self.slDispMin.valueChanged.connect(self.dispslider_changed)

        self.slDispMax.setMinimum(0)  # Set the minimum value
        self.slDispMax.setMaximum(int(self.max * 100))
        self.slDispMax.setSingleStep(int(self.max * 10))
        self.slDispMax.setValue(int(self.max * 100))
        self.lbDispMax.setText(f'{self.max: .2f}')

        self.slDispMin.setMinimum(int(self.min * 100))  # Set the minimum value
        self.slDispMin.setMaximum(0)
        self.slDispMin.setSingleStep(int(self.min * 10))
        self.slDispMin.setValue(int(self.min * 100))
        self.lbDispMin.setText(f'{self.min: .2f}')

        self.slfft = QSlider(self)
        self.slfft.setGeometry(680, 650, 291, 16)
        self.slfft.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.slfft.setObjectName("slfft")
        self.slfft.setEnabled(False)
        self.programmatic_change = True
        self.slfft.setMinimum(0)  # Set the minimum value
        self.slfft.setMaximum(1)
        self.slfft.setValue(1)
        self.programmatic_change = False
        self.slfft.valueChanged.connect(self.fftslider_changed)

        self.rect = None
        self.ROItext = None
        self.dispshow_image()
        self.select_roi()

    def dispshow_image(self):
        self.dispax.clear()
        min = float(self.slDispMin.value()) / 100
        max = float(self.slDispMax.value()) / 100
        extent = [self.imroixi, self.imroixf, self.imroizf, self.imroizi]
        image = self.dispax.imshow(self.dispmap[:, :, self.index - 1], cmap='jet',
                                   extent=extent, vmin = min, vmax = max, interpolation='none')


        self.dispax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispax.tick_params(axis='x', labelsize=9)
        self.dispax.tick_params(axis='y', labelsize=9)
        self.dispcanvas.draw()
        self.dispax.set_aspect('equal')
        self.generate_roi()

    def dispupdate_frame(self):
        min = float(self.slDispMin.value()) / 100
        max = float(self.slDispMax.value()) / 100
        self.lbDispMax.setText(f'{max: .2f}')
        self.lbDispMin.setText(f'{min: .2f}')
        self.dispshow_image()

    def dispslider_changed(self):
        self.dispupdate_frame()

    def select_roi(self):
        self.PSD = True
        if self.rect is not None:
            self.rect.remove()
            self.rect = None
            self.dispcanvas.draw()
        if self.ROItext is not None:
            self.ROItext.remove()
            self.dispcanvas.draw()
        self.ROItext = self.dispax.text(0.5, 0.5, "Select a ROI \n Right-Click to finish",
                                    color='black', fontsize=12, ha='center', va='center', transform=self.dispax.transAxes, )
        self.dispcanvas.draw()
        self.rs = RectangleSelector(self.dispax, self.draw_roi, useblit=True,
                                    button=[1], minspanx=5, minspany=5, spancoords='pixels',
                                    interactive=True,  props = dict(facecolor='black', edgecolor='black', alpha=0.3, fill=True))
        self.is_drawing = True
        self.connection_id = self.dispfig.canvas.mpl_connect('button_press_event', self.finish_roi)
        self.dispAxis.setCursor(Qt.CursorShape.CrossCursor)

    def draw_roi(self, eclick, erelease, *args, **kwargs):
        if self.is_drawing:
            if self.rect is not None:
                self.rect.remove()
                self.rect = None
                self.dispcanvas.draw()
            # Draw a red rectangle on the selected region
            self.rect = Rectangle((eclick.xdata, eclick.ydata),
                             erelease.xdata - eclick.xdata,
                             erelease.ydata - eclick.ydata,
                             linewidth=2, edgecolor='black', facecolor='none', linestyle = 'dashed')
            self.dispax.add_patch(self.rect)
            self.dispcanvas.draw()
            #

    def finish_roi(self,event):
        #if event.dblclick and event.button == 1:
        if event.button == 3:
            if self.rect is not None:
                self.rs.set_visible(False)
                self.rs.set_active(False)
                self.rs.disconnect_events()
                self.rs = None
                self.roixi = self.rect.get_x()
                self.roizi = self.rect.get_y()
                self.roixf = self.roixi + self.rect.get_width()
                self.roizf = self.roizi + self.rect.get_height()
                self.ROItext.remove()
                self.ROItext = None
                self.rect.remove()
                self.rect = None
                self.dispcanvas.draw()
                self.dispfig.canvas.mpl_disconnect(self.connection_id)
                self.is_drawing = False
                self.generate_roi()
                self.dispAxis.setCursor(Qt.CursorShape.ArrowCursor)
                self.dispshow_image()
                self.generate_signal()
                self.slDisp.setEnabled(True)
                self.slfft.setEnabled(True)
            else:
                return

    def generate_roi(self):
        self.rect = Rectangle((self.roixi, self.roizi),
        self.roixf - self.roixi, self.roizf - self.roizi,
                      linewidth=2, edgecolor='black', facecolor='none', linestyle='dashed')
        self.dispax.add_patch(self.rect)
        self.dispcanvas.draw()

    def generate_signal(self):
        self.signalMean = np.zeros((self.frames))
        self.timedim = np.zeros((self.frames))
        for f in range(0,self.frames):
            sum = 0
            count = 0
            for z in range(int((self.roizi-self.imroizi)/self.dz),int((self.roizf-self.imroizi)/self.dz)):
                for x in range(int((self.roixi-self.imroixi)/self.dx),int((self.roixf-self.imroixi)/self.dx)):
                    if ((0 <= z < self.Z) and (0 <= x < self.X)):
                        sum = sum + self.dispmap[z,x,f]
                        count = count + 1
            if count>0:
                self.signalMean[f] = sum/count
            self.timedim[f] = f*self.dt*1000
        self.plotSignal()

    def plotSignal(self):
        self.signalax.clear()
        self.signalax.set_facecolor((1, 1, 1))
        self.signalax.plot(self.timedim, self.signalMean, 'k-')
        self.signalcanvas.draw()

        self.plotMarker = self.signalax.plot(self.timedim[self.index - 1], self.signalMean[self.index - 1], 'ro')
        self.signalax.set_ylabel('Amplitude [μm]', fontsize=9)
        self.signalax.set_xlabel('Time [ms]', fontsize=9)
        self.signalax.tick_params(axis='x', labelsize=9)
        self.signalax.tick_params(axis='y', labelsize=9)
        self.signalax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        self.signalax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')
        self.signalax.set_title(f'{self.timedim[self.index - 1]: .3f} ms , {self.signalMean[self.index - 1]: .3f} μm')
        self.signalcanvas.draw()
        if self.PSD:
            self.plotPSD()

    def plotSignalupdate(self):
        self.signalax.lines[1].remove()
        self.plotMarker = self.signalax.plot(self.timedim[self.index - 1], self.signalMean[self.index - 1], 'ro')
        self.signalax.set_title(f'{self.timedim[self.index - 1]: .3f} ms , {self.signalMean[self.index - 1]: .3f} μm')
        self.signalcanvas.draw()

    def plotPSD(self):
        [self.ps, self.freq] = self.powerSpectrum(self.signalMean, 1/self.dt)
        self.psax.clear()
        self.psfig.set_facecolor((0.95, 0.95, 0.95))
        self.psax.set_facecolor((1, 1, 1))
        self.psax.plot(self.freq, self.ps, 'k-')
        self.psax.set_ylabel('Amplitude [a.u.]', fontsize=9)
        self.psax.set_xlabel('Frequency [Hz]', fontsize=9)
        self.psax.tick_params(axis='x', labelsize=9)
        self.psax.tick_params(axis='y', labelsize=9)
        self.psax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        self.psax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')
        self.fftplotMarker = self.psax.plot(self.freq[0], self.ps[0], 'ro')
        self.psax.set_title(
                f'{self.freq[0]: .3f} Hz , {self.ps[0]: .3f} a.u.')
        self.pscanvas.draw()
        self.PSD = False

    def update_frame(self):
        self.index = self.slDisp.value()
        self.lbDisp.setText(f'{self.index} / {self.frames} ')
        self.dispshow_image()
        #self.plotSignal()
        self.plotSignalupdate()

    def slider_changed(self):
        if self.programmatic_change:
            return
        else:
            self.update_frame()

    def fftslider_changed(self):
        if self.programmatic_change:
            return
        else:
            fftindex = self.slfft.value()
            self.psax.lines[1].remove()
            self.fftplotMarker = self.psax.plot(self.freq[fftindex], self.ps[fftindex], 'ro')
            self.psax.set_title(
                f'{self.freq[fftindex]: .3f} Hz , {self.ps[fftindex]: .3f} a.u.')
            self.pscanvas.draw()


    def powerSpectrum(self,signal, fs):
        # Compute the FFT of the signal
        interfft = np.fft.fft(signal, n=1024)

        # Compute the power spectrum (magnitude squared)
        interps = np.abs(interfft) ** 2

        # Frequencies corresponding to the FFT result
        # For real input, the FFT result is symmetric, and the positive frequencies are given by:
        frequencies = np.fft.fftfreq(len(interfft), d=1 / fs)
        positive_frequencies = frequencies[:len(frequencies) // 2]
        ps = interps[:len(interps) // 2]
        self.programmatic_change = True
        self.slfft.setMinimum(0)  # Set the minimum value
        self.slfft.setMaximum(positive_frequencies.shape[0]-1)
        self.programmatic_change = False
        return ps, positive_frequencies


    def closeEvent(self, event):
        return
       # self.retranslateUi(MainWindow)

class Velo_area(QDialog):

    def __init__(self, dispmap, prm, parent=None):
        super().__init__()
        self.initUI(dispmap,prm)

    def initUI(self, velomap, prm):
        self.velomap = velomap
        [self.Z, self.X] = (self.velomap.shape)

        self.dx = prm['dx']
        self.dz = prm['dz']

        self.roixi = 0
        self.roizi = 0
        self.roixf = self.X * self.dx
        self.roizf = self.Z * self.dz

        self.setWindowTitle("Velocity Area Analyzer")
        self.resize(730, 700)

        self.max = np.amax(self.velomap)
        self.min = np.amin(self.velomap)

        self.font = QFont()
        self.font.setPointSize(16)

        #Disp Axis
        self.veloAxis = QWidget(self)
        self.veloAxis.setObjectName("veloAxis")
        self.veloAxis.setGeometry(30, 10, 591, 611)
        layout = QVBoxLayout(self.veloAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.velofig = plt.Figure()
        self.velofig.set_facecolor((0.94, 0.94, 0.94))
        self.veloax = self.velofig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.veloax.set_xlabel('Lateral [mm]', fontsize=9)
        self.veloax.set_ylabel('Axial [mm]', fontsize=9)
        self.veloax.tick_params(axis='x', labelsize=9)
        self.veloax.tick_params(axis='y', labelsize=9)
        self.velocanvas = FigureCanvas(self.velofig)
        layout.addWidget(self.velocanvas)
        layout.addStretch()

        self.pbResetROI = QPushButton(self)
        self.pbResetROI.setGeometry(50, 640, 75, 23)
        self.pbResetROI.setObjectName("pbResetROI")
        self.pbResetROI.setText("Reset ROI")
        self.pbResetROI.clicked.connect(self.select_roi)

        self.lbMean = QLabel(self)
        self.lbMean.setGeometry(20, 10, 251, 40)
        self.lbMean.setObjectName("lbMean")
        self.lbMean.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbMean.setText(f'V<html><sub>ROI</sub></html> = -- m/s')
        self.lbMean.setFont(self.font)

        self.lbMax = QLabel(self)
        self.lbMax.setGeometry(610, 20, 81, 20)
        self.lbMax.setObjectName("lbMax")
        self.lbMax.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbMax.setText(f'{self.max : .2f} m/s')

        self.slMax = QSlider(self)
        self.slMax.setOrientation(Qt.Orientation.Vertical)
        self.slMax.setGeometry(600, 30, 16, 270)
        self.slMax.setObjectName("slMax")
        self.slMax.setMinimum(0)  # Set the minimum value
        self.slMax.setMaximum(int(self.max * 100))
        self.slMax.setSingleStep(int(10))
        self.slMax.setValue(int(self.max * 100))
        self.slMax.valueChanged.connect(self.veloslider_changed)

        self.lbMin = QLabel(self)
        self.lbMin.setGeometry(610, 550, 81, 20)
        self.lbMin.setObjectName("lbMax")
        self.lbMin.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbMin.setText(f'{self.min : .2f} m/s')

        self.slMin = QSlider(self)
        self.slMin.setOrientation(Qt.Orientation.Vertical)
        self.slMin.setGeometry(600, 300, 16, 270)
        self.slMin.setObjectName("slMin")
        self.slMin.setMinimum(-int(abs(self.min) * 100))  # Set the minimum value
        self.slMin.setMaximum(0)
        self.slMin.setSingleStep(int(10))
        self.slMin.setValue(-int(abs(self.min) * 100))
        #self.lbVeloMin.setText(f'{-self.gtr: .2f}')
        self.slMin.valueChanged.connect(self.veloslider_changed)

        self.rect = None
        self.ROItext = None
        self.veloshow_image()
        self.select_roi()

    def veloshow_image(self):
        self.veloax.clear()
        min = float(self.slMin.value()) / 100
        max = float(self.slMax.value()) / 100
        extent = [0, self.X*self.dx, self.Z*self.dz, 0]
        image = self.veloax.imshow(self.velomap, cmap='jet',
                                   extent=extent, vmin = min, vmax = max, interpolation='none')

        self.veloax.set_xlabel('Lateral [mm]', fontsize=9)
        self.veloax.set_ylabel('Axial [mm]', fontsize=9)
        self.veloax.tick_params(axis='x', labelsize=9)
        self.veloax.tick_params(axis='y', labelsize=9)
        self.velocanvas.draw()
        self.veloax.set_aspect('equal')
        self.generate_roi()

    def veloupdate_frame(self):
        min = float(self.slMin.value()) / 100
        max = float(self.slMax.value()) / 100
        self.lbMax.setText(f'{max : .2f} m/s')
        self.lbMin.setText(f'{min : .2f} m/s')

        #self.lbVeloMax.setText(f'{max: .2f}')
        #self.lbVeloMin.setText(f'{min: .2f}')
        self.veloshow_image()

    def veloslider_changed(self):
        self.veloupdate_frame()

    def select_roi(self):
        if self.rect is not None:
            self.rect.remove()
            self.rect = None
            self.lbMean.setText(f'V<html><sub>ROI</sub></html> = -- m/s')
            self.velocanvas.draw()
        if self.ROItext is not None:
            self.ROItext.remove()
            self.velocanvas.draw()
        self.ROItext = self.veloax.text(0.5, 0.5, "Select a ROI \n Right-Click to finish",
                                    color='black', fontsize=12, ha='center', va='center', transform=self.veloax.transAxes, )
        self.velocanvas.draw()
        self.rs = RectangleSelector(self.veloax, self.draw_roi, useblit=True,
                                    button=[1], minspanx=5, minspany=5, spancoords='pixels',
                                    interactive=True,  props = dict(facecolor='black', edgecolor='black', alpha=0.3, fill=True))
        self.is_drawing = True
        self.connection_id = self.velofig.canvas.mpl_connect('button_press_event', self.finish_roi)
        self.veloAxis.setCursor(Qt.CursorShape.CrossCursor)

    def draw_roi(self, eclick, erelease, *args, **kwargs):
        if self.is_drawing:
            if self.rect is not None:
                self.rect.remove()
                self.rect = None
                self.velocanvas.draw()
            # Draw a red rectangle on the selected region
            self.rect = Rectangle((eclick.xdata, eclick.ydata),
                             erelease.xdata - eclick.xdata,
                             erelease.ydata - eclick.ydata,
                             linewidth=2, edgecolor='black', facecolor='none', linestyle = 'dashed')
            self.veloax.add_patch(self.rect)
            self.velocanvas.draw()
            #

    def finish_roi(self,event):
        #if event.dblclick and event.button == 1:
        if event.button == 3:
            if self.rect is not None:
                self.rs.set_visible(False)
                self.rs.set_active(False)
                self.rs.disconnect_events()
                self.rs = None
                self.roixi = self.rect.get_x()
                self.roizi = self.rect.get_y()
                self.roixf = self.roixi + self.rect.get_width()
                self.roizf = self.roizi + self.rect.get_height()
                self.ROItext.remove()
                self.ROItext = None
                self.rect.remove()
                self.rect = None
                self.velocanvas.draw()
                self.velofig.canvas.mpl_disconnect(self.connection_id)
                self.is_drawing = False
                self.generate_roi()
                self.veloAxis.setCursor(Qt.CursorShape.ArrowCursor)
                self.veloshow_image()
                self.generate_signal()
            else:
                return

    def generate_roi(self):
        self.rect = Rectangle((self.roixi, self.roizi),
        self.roixf - self.roixi, self.roizf - self.roizi,
                      linewidth=2, edgecolor='black', facecolor='none', linestyle='dashed')
        self.veloax.add_patch(self.rect)
        self.velocanvas.draw()

    def generate_signal(self):
        self.signalMean = 0
        self.signalStd = 0
        count = 1
        for z in range(int(self.roizi/self.dz),int(self.roizf/self.dz)):
            for x in range(int(self.roixi/self.dx),int(self.roixf/self.dx)):
                if((0 <= z < self.Z) and (0 <= x < self.X)):
                    if self.velomap[z,x] != 0:
                        self.signalMean = self.signalMean + self.velomap[z,x]
                        count = count + 1

        self.signalMean = self.signalMean/count

        # Calculate standard deviation
        sum_squared_diff = 0

        for z in range(int(self.roizi / self.dz), int(self.roizf / self.dz)):
            for x in range(int(self.roixi / self.dx), int(self.roixf / self.dx)):
                if self.velomap[z, x] != 0:
                    if ((0 <= z < self.Z) and (0 <= x < self.X)):
                        diff = self.velomap[z, x] - self.signalMean
                        sum_squared_diff += diff ** 2

        self.signalStd = np.sqrt(sum_squared_diff / count)

        self.lbMean.setText(f' V<html><sub>ROI</sub></html> = {self.signalMean: .2f} ±  {self.signalStd: .2f} m/s')

    def closeEvent(self, event):
        return
       # self.retranslateUi(MainWindow)

class phase_velo(QDialog):
    dispcurvedata = pyqtSignal(numpy.ndarray,numpy.ndarray)
    def __init__(self, dispmap, prm, parent=None):
        super().__init__()
        self.initUI(dispmap,prm)

    def initUI(self, dispmap, prm):
        self.dispmap = dispmap
        [self.Z, self.X, self.frames] = (self.dispmap.shape)

        self.dx = prm['dx']
        self.dz = prm['dz']
        self.dt = prm['dt']
        self.index = prm['index']
        self.zeropad = prm['zeropad']
        aux = prm['tr']
        self.tr = pow(10,(float(aux)/20))
        self.freqcut = prm['freqcut']
        self.df = prm['df']
        self.density = prm['density']
        self.fs = 1/self.dt

        self.roixi = 0
        self.roizi = 0
        self.roixf = self.X * self.dx
        self.roizf = self.Z * self.dz

        self.setWindowTitle("Phase Velocity Analyzer")
        self.resize(1000, 700)

        self.max = np.amax(self.dispmap)
        self.min = np.amin(self.dispmap)

        self.font = QFont()
        self.font.setPointSize(16)

        #Disp Axis
        self.dispmapAxis = QWidget(self)
        self.dispmapAxis.setObjectName("dispmapAxis")
        self.dispmapAxis.setGeometry(20, 20, 391, 311)
        layout = QVBoxLayout(self.dispmapAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.dispfig = plt.Figure()
        self.dispfig.set_facecolor((0.94, 0.94, 0.94))
        self.dispax = self.dispfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.dispax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispax.tick_params(axis='x', labelsize=9)
        self.dispax.tick_params(axis='y', labelsize=9)
        self.dispcanvas = FigureCanvas(self.dispfig)
        layout.addWidget(self.dispcanvas)
        layout.addStretch()

        # Disp Slider
        self.slDisp = QSlider(self)
        self.slDisp.setGeometry(30, 340, 371, 16)
        self.slDisp.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.slDisp.setObjectName("slDisp")
        self.slDisp.setEnabled(True)
        self.programmatic_change = True
        self.slDisp.setMinimum(1)  # Set the minimum value
        self.slDisp.setMaximum(self.frames)
        self.slDisp.setValue(self.index)
        self.programmatic_change = False
        self.slDisp.valueChanged.connect(self.slider_changed)

        self.slMax = QSlider(self)
        self.slMax.setGeometry(400, 20, 16, 151)
        self.slMax.setOrientation(QtCore.Qt.Orientation.Vertical)
        self.slMax.setObjectName("slMax")
        self.slMax.setEnabled(True)
        self.programmatic_change = True
        self.slMax.setMinimum(0)  # Set the minimum value
        self.slMax.setMaximum(int(self.max*100))
        self.slMax.setValue(int(self.max*100))
        self.programmatic_change = False
        self.slMax.valueChanged.connect(self.slider_changed)

        self.slMin = QSlider(self)
        self.slMin.setGeometry(400, 170, 16, 151)
        self.slMin.setOrientation(QtCore.Qt.Orientation.Vertical)
        self.slMin.setObjectName("slMin")
        self.slMin.setEnabled(True)
        self.programmatic_change = True
        self.slMin.setMinimum(int(self.min * 100))  # Set the minimum value
        self.slMin.setMaximum(0)
        self.slMin.setValue(int(self.min * 100))
        self.programmatic_change = False
        self.slMin.valueChanged.connect(self.slider_changed)

        self.lbMax = QLabel(self)
        self.lbMax.setGeometry(430, 20, 61, 20)
        self.lbMax.setObjectName("lbMax")
        self.lbMin = QLabel(self)
        self.lbMin.setGeometry(430, 310, 61, 20)
        self.lbMin.setObjectName("lbMin")


        # Proj Axis
        self.projAxis = QWidget(self)
        self.projAxis.setObjectName("projAxis")
        self.projAxis.setGeometry(490, 10, 391, 311)
        layout = QVBoxLayout(self.projAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.projfig = plt.Figure()
        self.projfig.set_facecolor((0, 0, 0))
        self.projax = self.projfig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.projax.set_facecolor((0, 0, 0))
        self.projax.set_xlabel('Lateral [mm]', fontsize=9)
        self.projax.set_ylabel('Axial [mm]', fontsize=9)
        self.projax.tick_params(axis='x', labelsize=9)
        self.projax.tick_params(axis='y', labelsize=9)
        self.projcanvas = FigureCanvas(self.projfig)
        layout.addWidget(self.projcanvas)
        layout.addStretch()

        # FFT Axis
        self.fftAxis = QWidget(self)
        self.fftAxis.setObjectName("fftAxis")
        self.fftAxis.setGeometry(20, 370, 391, 311)
        layout = QVBoxLayout(self.fftAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.fftfig = plt.Figure()
        self.fftfig.set_facecolor((0, 0, 0))
        self.fftax = self.fftfig.add_axes(rect=[0.20, 0.14, 0.75, 0.8])
        self.fftax.set_facecolor((0, 0, 0))
        self.fftax.set_xlabel('Lateral [mm]', fontsize=9)
        self.fftax.set_ylabel('Axial [mm]', fontsize=9)
        self.fftax.tick_params(axis='x', labelsize=9)
        self.fftax.tick_params(axis='y', labelsize=9)
        self.fftcanvas = FigureCanvas(self.fftfig)
        layout.addWidget(self.fftcanvas)
        layout.addStretch()

        # dispcurve Axis
        self.dispcurveAxis = QWidget(self)
        self.dispcurveAxis.setObjectName("dispcurveAxis")
        self.dispcurveAxis.setGeometry(490, 370, 391, 311)
        layout = QVBoxLayout(self.dispcurveAxis)
        layout.setContentsMargins(1, 1, 1, 1)
        layout.addStretch()
        self.dispcurvefig = plt.Figure()
        self.dispcurvefig.set_facecolor((0, 0, 0))
        self.dispcurveax = self.dispcurvefig.add_axes(rect=[0.17, 0.14, 0.75, 0.8])
        self.dispcurveax.set_facecolor((0, 0, 0))
        self.dispcurveax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispcurveax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispcurveax.tick_params(axis='x', labelsize=9)
        self.dispcurveax.tick_params(axis='y', labelsize=9)
        self.dispcurvecanvas = FigureCanvas(self.dispcurvefig)
        layout.addWidget(self.dispcurvecanvas)
        layout.addStretch()

        self.pbReset = QPushButton(self)
        self.pbReset.setGeometry(900, 40, 75, 23)
        self.pbReset.setObjectName("pbReset")
        self.pbReset.setText("Reset")
        self.pbReset.clicked.connect(self.select_roi)

        self.pbDone = QPushButton(self)
        self.pbDone.setGeometry(900, 80, 75, 23)
        self.pbDone.setObjectName("pbDone")
        self.pbDone.setText("Done")
        self.pbDone.setEnabled(False)
        self.pbDone.clicked.connect(self.finish)

        self.rect = None
        self.ROItext = None
        self.phrect = None
        self.phROItext = None
        self.projected = False
        self.projROI =True
        self.update_frame()
        self.select_roi()

    def dispshow_image(self):
        self.dispax.clear()
        min = float(self.slMin.value()) / 100
        max = float(self.slMax.value()) / 100
        extent = [0, self.X*self.dx, self.Z*self.dz, 0]
        #image = self.dispax.imshow(self.dispmap[:,:,self.index-1], cmap='jet',
        #                           extent=extent, vmin = min, vmax = max, interpolation='none', aspect='auto')
        image = self.dispax.imshow(self.dispmap[:, :, self.index - 1], cmap='jet',
                                    extent=extent, vmin = min, vmax = max, interpolation='none')

        self.dispax.set_xlabel('Lateral [mm]', fontsize=9)
        self.dispax.set_ylabel('Axial [mm]', fontsize=9)
        self.dispax.tick_params(axis='x', labelsize=9)
        self.dispax.tick_params(axis='y', labelsize=9)
        self.dispax.set_title(f'Frame # {self.index : .0f} ', fontsize=9)
        self.dispcanvas.draw()
        self.generate_roi()

    def proj_show(self):
        if self.projected:
            self.projfig.set_facecolor((0.94, 0.94, 0.94))
            self.projax.clear()
            min = float(self.slMin.value()) / 100
            max = float(self.slMax.value()) / 100
            extent = [0, self.frames*self.dt*1000, self.X*self.dx,0 ]
            image = self.projax.imshow(self.proj, cmap='jet',
                                   extent=extent, vmin = min, vmax = max, interpolation='none', aspect='auto')

            self.projax.set_xlabel('Time [ms]', fontsize=9)
            self.projax.set_ylabel('Lateral [mm]', fontsize=9)
            self.projax.tick_params(axis='x', labelsize=9)
            self.projax.tick_params(axis='y', labelsize=9)

            self.projcanvas.draw()
            if self.projROI:
                self.phselect_roi()
        else:
            return

    def fft_show(self):
        if self.projected:
            self.fftfig.set_facecolor((0.94, 0.94, 0.94))
            self.fftax.clear()
            extent = [self.dimf[0], self.dimf[-1], self.dimkx[-1],self.dimkx[0]]
            image = self.fftax.imshow(np.abs(self.fftdata), cmap='inferno',
                                   extent=extent, interpolation='none', aspect='auto')
            data = np.abs(self.fftdata)
            #image = self.fftax.imshow(data , cmap='inferno',
            #                          extent=extent, interpolation='none', aspect='auto',norm=mcolors.PowerNorm(gamma=1))

            self.fftax.set_xlabel('Frequency [Hz]', fontsize=9)
            self.fftax.set_ylabel('k [mm-1]', fontsize=9)
            self.fftax.tick_params(axis='x', labelsize=9)
            self.fftax.tick_params(axis='y', labelsize=9)
            self.fftax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            self.fftcanvas.draw()
        else:
            return

    def dispcurve_show(self):
        if self.projected:
            self.dispcurvefig.set_facecolor((0.94, 0.94, 0.94))
            self.dispcurveax.clear()
            self.dispcurveax.set_facecolor((1, 1, 1))
            self.dispcurveax.plot(self.dimf, self.velo, marker='None', linestyle='-', color='k')

            self.dispcurveax.set_xlabel('Frequency [Hz]', fontsize=9)
            self.dispcurveax.set_ylabel('Velocity [m/s]', fontsize=9)
            self.dispcurveax.tick_params(axis='x', labelsize=9)
            self.dispcurveax.tick_params(axis='y', labelsize=9)
            self.dispcurveax.grid(which='both', linestyle='-', linewidth='0.5', color='gray')
            self.dispcurveax.set_xlim([0, self.freqcut])
            self.dispcurveax.set_ylim([0, 6])
            self.dispcurvecanvas.draw()
            self.fftax.plot(self.dimf, self.kmax, marker='o',color='g',markersize=0.5)
            self.fftax.set_xlim([-self.freqcut, self.freqcut])
            self.fftax.set_ylim([-self.freqcut, self.freqcut])
            self.fftax.figure.canvas.draw()
            self.pbDone.setEnabled(True)
        else:
            return

    def update_frame(self):
        self.index = self.slDisp.value()
        min = float(self.slMin.value()) / 100
        max = float(self.slMax.value()) / 100
        self.lbMax.setText(f'{max : .2f} ')
        self.lbMin.setText(f'{min : .2f} ')

        #self.lbVeloMax.setText(f'{max: .2f}')
        #self.lbVeloMin.setText(f'{min: .2f}')
        self.dispshow_image()
        self.proj_show()

    def slider_changed(self):
        self.update_frame()

    def select_roi(self):
        if self.rect is not None:
            self.rect.remove()
            self.rect = None
            self.dispcanvas.draw()
        if self.ROItext is not None:
            self.ROItext.remove()
            self.dispcanvas.draw()
        self.ROItext = self.dispax.text(0.5, 0.5, "Select a ROI \n Right-Click to finish",
                                    color='black', fontsize=12, ha='center', va='center', transform=self.dispax.transAxes, )
        self.dispcanvas.draw()
        self.rs = RectangleSelector(self.dispax, self.draw_roi, useblit=True,
                                    button=[1], minspanx=5, minspany=5, spancoords='pixels',
                                    interactive=True,  props = dict(facecolor='black', edgecolor='black', alpha=0.3, fill=True))
        self.is_drawing = True
        self.projROI = True
        self.connection_id = self.dispfig.canvas.mpl_connect('button_press_event', self.finish_roi)
        self.dispmapAxis.setCursor(Qt.CursorShape.CrossCursor)

    def draw_roi(self, eclick, erelease, *args, **kwargs):
        if self.is_drawing:
            if self.rect is not None:
                self.rect.remove()
                self.rect = None
                self.dispcanvas.draw()
            # Draw a red rectangle on the selected region
            self.rect = Rectangle((eclick.xdata, eclick.ydata),
                             erelease.xdata - eclick.xdata,
                             erelease.ydata - eclick.ydata,
                             linewidth=2, edgecolor='black', facecolor='none', linestyle = 'dashed')
            self.dispax.add_patch(self.rect)
            self.dispcanvas.draw()
            #

    def finish_roi(self,event):
        #if event.dblclick and event.button == 1:
        if event.button == 3:
            if self.rect is not None:
                self.rs.set_visible(False)
                self.rs.set_active(False)
                self.rs.disconnect_events()
                self.rs = None
                self.roixi = self.rect.get_x()
                self.roizi = self.rect.get_y()
                self.roixf = self.roixi + self.rect.get_width()
                self.roizf = self.roizi + self.rect.get_height()
                self.ROItext.remove()
                self.ROItext = None
                self.rect.remove()
                self.rect = None
                self.dispcanvas.draw()
                self.dispfig.canvas.mpl_disconnect(self.connection_id)
                self.is_drawing = False
                self.generate_roi()
                self.dispmapAxis.setCursor(Qt.CursorShape.ArrowCursor)
                self.dispshow_image()
                self.generate_signal()
            else:
                return

    def generate_roi(self):
        self.rect = Rectangle((self.roixi, self.roizi),
        self.roixf - self.roixi, self.roizf - self.roizi,
                      linewidth=2, edgecolor='black', facecolor='none', linestyle='dashed')
        self.dispax.add_patch(self.rect)
        self.dispcanvas.draw()

    def generate_signal(self):
        self.proj = np.zeros([self.X,self.frames])
        roixlim1 = (self.roixi/self.dx)
        roixlim2 = (self.roixf / self.dx)
        roizlim1 = (self.roizi / self.dz)
        roizlim2 = (self.roizf / self.dz)

        for frame in range(0, self.frames):
            for x in range(0,self.X):
                if (roixlim1 <= x < roixlim2):
                    count = 0
                    sum1 = 0
                    for z in range(0,self.Z):
                        if (roizlim1<=z<roizlim2):
                            sum1 = sum1 + self.dispmap[z,x,frame]
                            count = count + 1
                    self.proj[x, frame] = sum1/count
        self.projected = True
        self.proj_show()

    def finish(self):
        lim1 = np.argmin(np.abs(self.dimf - 0))
        lim2 = np.argmin(np.abs(self.dimf - self.freqcut))

        index = 0
        dispdata = np.zeros([lim2 - lim1])
        freqdata = np.zeros([lim2 - lim1])
        for i in range(lim1,lim2):
            dispdata[index] = self.velo[i]
            freqdata[index] = self.dimf[i]
            index = index+1

        self.dispcurvedata.emit(freqdata,dispdata)

        self.accept()

    def phselect_roi(self):
        if self.phrect is not None:
            self.phrect.remove()
            self.phrect = None
            self.projcanvas.draw()
        if self.phROItext is not None:
            self.phROItext.remove()
            self.projcanvas.draw()
        self.phROItext = self.projax.text(0.5, 0.5, "Select a ROI \n Right-Click to finish",
                                    color='black', fontsize=12, ha='center', va='center', transform=self.projax.transAxes, )
        self.projcanvas.draw()
        self.phrs = RectangleSelector(self.projax, self.phdraw_roi, useblit=True,
                                    button=[1], minspanx=5, minspany=5, spancoords='pixels',
                                    interactive=True,  props = dict(facecolor='black', edgecolor='black', alpha=0.3, fill=True))
        self.phis_drawing = True
        self.phconnection_id = self.projfig.canvas.mpl_connect('button_press_event', self.phfinish_roi)
        self.projAxis.setCursor(Qt.CursorShape.CrossCursor)

    def phdraw_roi(self, eclick, erelease, *args, **kwargs):
        if self.phis_drawing:
            if self.phrect is not None:
                self.phrect.remove()
                self.phrect = None
                self.projcanvas.draw()
            # Draw a red rectangle on the selected region
            self.phrect = Rectangle((eclick.xdata, eclick.ydata),
                             erelease.xdata - eclick.xdata,
                             erelease.ydata - eclick.ydata,
                             linewidth=2, edgecolor='black', facecolor='none', linestyle = 'dashed')
            self.projax.add_patch(self.phrect)
            self.projcanvas.draw()
            #

    def phfinish_roi(self,event):
        #if event.dblclick and event.button == 1:
        if event.button == 3:
            if self.phrect is not None:
                self.phrs.set_visible(False)
                self.phrs.set_active(False)
                self.phrs.disconnect_events()
                self.phrs = None
                self.phroixi = self.phrect.get_x()
                self.phroizi = self.phrect.get_y()
                self.phroixf = self.phroixi + self.phrect.get_width()
                self.phroizf = self.phroizi + self.phrect.get_height()
                self.phROItext.remove()
                self.phROItext = None
                self.phrect.remove()
                self.phrect = None
                self.projcanvas.draw()
                self.projfig.canvas.mpl_disconnect(self.phconnection_id)
                self.phis_drawing = False
                self.phgenerate_roi()
                self.projAxis.setCursor(Qt.CursorShape.ArrowCursor)
                self.phgenerate_roi()
                self.phgenerate_image()
                self.projROI = False

            else:
                return

    def phgenerate_roi(self):
        self.phrect = Rectangle((self.phroixi, self.phroizi),
        self.phroixf - self.phroixi, self.phroizf - self.phroizi,
                      linewidth=2, edgecolor='black', facecolor='none', linestyle='dashed')
        self.projax.add_patch(self.phrect)
        self.projcanvas.draw()

    def phgenerate_image(self):

        roixlim1 = (self.phroixi / self.dt)/1000
        roixlim2 = (self.phroixf / self.dt)/1000
        roizlim1 = (self.phroizi / self.dx)
        roizlim2 = (self.phroizf / self.dx)
        projmask = np.zeros([int(abs(roizlim2-roizlim1))+1,int(abs(roixlim2-roixlim1))+1])

        cx = 0
        for x in range(0, self.frames):
            cz = 0
            if (roixlim1 <= x < roixlim2):
                for z in range(0,self.X):
                    if (roizlim1 <= z < roizlim2):
                        projmask[cz, cx] = self.proj[z,x]
                        cz = cz + 1
                cx = cx + 1
        #ifft_1 = np.fft.ifft(np.flipup(projmask), axis=0, n= 2048)
        ifft_1 = np.fft.ifft(projmask, axis=0, n= self.zeropad)
        ifft_2 = np.fft.ifft(ifft_1, axis=1, n= self.zeropad)

        # Calculate df
        fkxdata_size = ifft_2.shape[0]
        df = self.fs / fkxdata_size

        # Calculate dimf
        self.dimf = np.linspace(-(self.fs / 2), (self.fs / 2), self.zeropad)

        # Calculate dimkx
        self.dimkx = np.linspace(-1 / (2 * (self.dx / 1000)), 1 / (2 * (self.dx / 1000)), self.zeropad)
        self.fftdata = ifft_2
        self.fftdata = np.fft.fftshift(self.fftdata, axes=0)
        self.fftdata = np.fft.fftshift(self.fftdata, axes=1)
        self.fftdata = np.flipud(self.fftdata)
        self.fft_show()
        self.findfftmax()

    def findfftmax(self):
        [Z,X] = self.fftdata.shape
        auxarray = np.zeros([Z])
        self.velo = np.zeros([X])
        self.kmax = np.zeros([X])
        fftmax = np.max(np.abs(self.fftdata))
        for x in range(0, X):
            auxarray = np.abs(self.fftdata[:,x])/fftmax
            if np.max(auxarray) > self.tr:
                for k in range(0,5):
                    p = np.argmax(auxarray)
                    self.velo[x] = self.velo[x]+abs(self.dimf[x]/self.dimkx[p])
                    self.kmax[x] = self.kmax[x]+self.dimkx[p]
                    auxarray[p] = 0
                self.velo[x] = self.velo[x] / 5
                self.kmax[x] = self.kmax[x] / 5
            else:
                self.velo[x] = 0
                self.kmax[x] = 0

        self.dispcurve_show()


    def closeEvent(self, event):
        self.dispcurvedata.emit(np.zeros([1]),np.zeros([1]))
       # self.retranslateUi(MainWindow)

class Error_Dialog(QDialog):

    def __init__(self, title):
        super().__init__()
        self.setWindowTitle(title)
        self.resize(300, 120)  # Adjust size
        # Create the label with an initial empty text
        self.label = QLabel('', self)
        self.label.setGeometry(20, 30, 260, 60)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        self.label.setWordWrap(True)
        # Create the "Finish" button
        self.pbFinish = QPushButton("OK", self)
        self.pbFinish.setGeometry(115, 80, 75, 30)
        self.pbFinish.clicked.connect(self.finish_clicked)

    def set_text(self, text):
        self.label.setText(text)

    def finish_clicked(self):
        self.accept()  # Close the dialog


    '''start_time = time.time()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")'''