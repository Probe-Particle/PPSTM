from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

import re
import sys

import os
import matplotlib
import numpy as np
matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from pyPPSTM.guiMethods import newPPSTM_simple, conv2float, importData
import pyPPSTM.basUtils as Bu
import pyPPSTM.elements as elements

######################### Canvas classes for plotting ####################

class SampleCanvas(FigureCanvasQTAgg):

    def __init__(self, parentWiget=None, parentApp=None,  width=5, height=4, dpi=100 ):
        self.fig  = Figure( figsize=(width, height), dpi=dpi )
        self.axes = self.fig.add_subplot(111)
                #super(self.__class__, self).__init__( self.fig )
        FigureCanvasQTAgg.__init__(self, self.fig )
        self.parent = parentApp
        self.setParent(parentWiget)
        FigureCanvasQTAgg.setSizePolicy(self,QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)
    
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parentApp=None, parentWidget=None, width=20, height=8, dpi=150):
        self.fig = plt.gcf()
        self.axes = self.fig.gca()
        super(MplCanvas, self).__init__(self.fig)
        self.parent = parentApp
        self.setParent(parentWidget)

    def correct_ext(self, fname, ext):
        _, fext = os.path.splitext(fname)
        if fext.capitalize() != ext.capitalize():
            fname += ext
        return fname
    
    def saveData(self):
        mapType = self.parent.map
        vv = self.parent.Vindx
        k = self.parent.Hindx
        namez = self.parent.plotData['namez']
        tip_type = self.parent.plotData['tip_type']
        WorkFunction = self.parent.plotData['WorkFunction']
        Voltages = self.parent.plotData['Voltages']
        WF_decay = self.parent.plotData['WF_decay']
        eta = self.parent.myDict['etaValue']
        if mapType == 'dIdV':
            fileName, fext = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()",'didv_'+namez[vv]+"_tip_"+tip_type+"-"+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d' %k, "WSxM Files (*.xyz);;XSF Files (*.xsf);;NPY Files (*.npy)")
        else:
            fileName, fext = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()",'STM_'+namez[vv]+"_tip_"+tip_type+"-"+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'_%03d' %k, "WSxM Files (*.xyz);;XSF Files (*.xsf);;NPY Files (*.npy)")

        if fileName:
            if fext == '(*.xsf)':
                print("For XSF or NPY outputs or tip_type = relaxed you have to have installed PPAFM in your PPSTM directory ")
                import pyPPSTM.GridUtils as GU
                print("writing XSF files")
                geom_plot = self.parent.plotData['geom_plot']
                lvec = self.parent.plotData['lvec']
                xsf_head = Bu.At2XSF(geom_plot) #if plot_atoms else GU.XSF_HEAD_DEFAULT
                if mapType == 'dIdV':
                    didv = self.parent.plotData['didv']
                    GU.saveXSF(fileName, didv[vv], lvec, head=xsf_head )
                elif mapType == 'STM':
                    current = self.parent.plotData['current']
                    GU.saveXSF(fileName, current[vv], lvec, head=xsf_head )
                print("XSF files written")

            elif fext == '(*.npy)':
                print("For XSF or NPY outputs or tip_type = relaxed you have to have installed PPAFM in your PPSTM directory ")
                import pyPPSTM.GridUtils as GU
                print("writing npy binary files")
                lvec = self.parent.plotData['lvec']
                if mapType == 'dIdV':
                    didv = self.parent.plotData['didv']
                    GU.saveNpy(fileName, didv[vv], lvec)#, head=XSF_HEAD_DEFAULT )
                elif mapType == 'STM':
                    current = self.parent.plotData['current']
                    GU.saveNpy(fileName, current[vv], lvec)#, head=XSF_HEAD_DEFAULT )
                print("npy files written")

            else:
                # default write in WSxM
                print("writing WSxM files")
                tip_r0 = self.parent.plotData['tip_r0']
                if mapType == 'dIdV':
                    didv = self.parent.plotData['didv']
                    tmp_curr=didv[vv,k,:,:].flatten()
                    out_curr=np.zeros((len(tmp_curr),3))
                    out_curr[:,0]=tip_r0[k,:,:,0].flatten()
                    out_curr[:,1]=tip_r0[k,:,:,1].flatten()
                    out_curr[:,2]=tmp_curr.copy()
                    f=open(fileName,'w')
                    print("WSxM file copyright Nanotec Electronica", file=f)
                    print("WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al.", file=f)
                    print("X[A]  Y[A]  Z[A]", file=f)
                    print("", file=f)
                    np.savetxt(f, out_curr)
                    f.close()
                elif mapType == 'STM':
                    current = self.parent.plotData['current']
                    tmp_curr=current[vv,k,:,:].flatten()
                    out_curr=np.zeros((len(tmp_curr),3))
                    out_curr[:,0]=tip_r0[k,:,:,0].flatten()
                    out_curr[:,1]=tip_r0[k,:,:,1].flatten()
                    out_curr[:,2]=tmp_curr.copy()
                    f=open(fileName,'w')
                    print("WSxM file copyright Nanotec Electronica", file=f)
                    print("WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al.", file=f)
                    print("X[A]  Y[A]  Z[A]", file=f)
                    print("", file=f)
                    np.savetxt(f, out_curr)
                    f.close()
                print("WSxM files written")

        print ("Data saved")

    def saveFigure(self):
        plotData = self.parent.plotData
        vv = self.parent.Vindx
        eta = self.parent.myDict['etaValue']
        k = self.parent.Hindx
        mapType = self.parent.map
        if mapType == 'dIdV':
            fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()",'didv_'+plotData['namez'][vv]+"_tip_"+plotData['tip_type']+"-"+"_WF_"+str(plotData['WorkFunction']-plotData['Voltages'][vv]*plotData['WF_decay'])+"_eta_"+str(eta)+'_%03d.png' %k,"Image files (*.png)")
            if fileName:
                fileName = self.correct_ext(fileName, ".png")
                print(("saving image to :", fileName))
                plt.fig = plt.gcf()
                plt.fig.savefig(fileName)
        else:
            fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()",'STM_'+plotData['namez'][vv]+"_tip_"+plotData['tip_type']+"-"+"_WF_"+str(plotData['WorkFunction'])+"_WF_decay_"+str(round(plotData['WF_decay'],1))+"_eta_"+str(eta)+'_%03d.png' %k,"Image files (*.png)")
            if fileName:
                fileName = self.correct_ext(fileName, ".png")
                print(("saving image to :", fileName))
                plt.savefig(fileName)

######################## Application Window Class ########################

class Window(QMainWindow):

    def __init__(self):
        super(Window, self).__init__()

        # Following are the variables from inputs. Kind of like private
        # Putting all the self variables to a dictionary
        # List of needed parameters:
        self.myDict = {'dft_code': 'fireball',
                       'sample_orbs': 'sp',
                       'spin': None,
                       'pbc': '00',
                       'data_format': 'xsf',
                       'kValue': 0.24,
                       'qValue': 0.0,
                       'x': [0.0,20.0,0.1],
                       'y': [0.0,20.0,0.1],
                       'z': [5.0,6.0,1.0],
                       'scan_type': 'dIdV',
                       'etaValue': 0.1,
                       'wf_decay': 0.0,
                       'V': -2.0,
                       'Vmax': 2.0,
                       'dV': 0.1,
                       'tipOrbS': 0.13,
                       'tipOrbPxy': 0.87,
                       'OMP_NUM_THREADS': 1,
                       'tip_type': 'fixed',}
            
        self.importData = None
        self.plotData = None

        self.string_widgets = {}
        self.num_widgets = {}
        self.paramList = ['dft_code', 'sample_orbs', 'spin', 'pbc', 'data_format', 'kValue', 'qValue', 'x', 'y', 'z', 'scan_type', 'etaValue', 'wf_decay', 'V', 'Vmax', 'dV', 'tipOrbS', 'tipOrbPxy', 'OMP_NUM_THREADS', 'tip_type', 'paths']
        self.runClicked = False
        self.paths = None
        self.map = 'dIdV'
        self.Hindx = 0
        self.Vindx = 0

        # App settings
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setWindowTitle("Application Window")

        # 'Dummy' widget to hold the whole layout
        self.mainWidget = QWidget(self)

        # Main layout. Horizontal box. Graphics on the left, functionality on the right.
        self.layout = QHBoxLayout(self.mainWidget)

        # Create a dummy widget for the plot
        canvLayout = QVBoxLayout()
        self.layout.addLayout(canvLayout)
        self.myCanvas = SampleCanvas(parentWiget=self.mainWidget, parentApp=self, width=5, height=4, dpi=100)
        canvLayout.addWidget(self.myCanvas)

        # Control Part on the right. Vertical Box Layout.
        controlLayout = QVBoxLayout()
        self.layout.addLayout(controlLayout)
        
####################################  INPUT PART  ##############################################

        # Input part 1 - HBox. First top region of the control part. 
        inputLayout1 = QHBoxLayout()
        controlLayout.addLayout(inputLayout1)

        # Adding codes to layout
        codesBox = QVBoxLayout(); inputLayout1.addLayout(codesBox)
        codesBox.addWidget(QLabel("codes"))

        # codes scroll bar
        self.codes = QComboBox(); codesBox.addWidget(self.codes)
        self.codes.addItems(['fireball', 'gpaw', 'aims', 'cp2k'])
        self.codes.currentIndexChanged[str].connect(self.selectCode)
        self.string_widgets['dft_code'] = self.codes
        
        # Adding orbitals to layout
        orbitalBox = QVBoxLayout(); inputLayout1.addLayout(orbitalBox)
        orbitalBox.addWidget(QLabel("sample orbital"))

        # oribital scroll bar
        self.orbitals = QComboBox(); orbitalBox.addWidget(self.orbitals)
        self.orbitals.addItems(['sp', 'spd'])
        self.orbitals.currentIndexChanged[str].connect(self.selectOrbital)
        self.string_widgets['sample_orbs'] = self.orbitals

        # Adding spin to Layout
        spinBox = QVBoxLayout(); inputLayout1.addLayout(spinBox)
        spinBox.addWidget(QLabel("spin"))

        # Spin scroll bar
        self.spin = QComboBox(); spinBox.addWidget(self.spin)
        self.spin.addItems(['None', 'both', 'alpha', 'beta'])
        self.spin.currentIndexChanged[str].connect(self.selectSpin)
        self.string_widgets['spin'] = self.spin

        # Adding pbc to layout
        pbcBox = QVBoxLayout(); inputLayout1.addLayout(pbcBox)
        pbcBox.addWidget(QLabel("pbc"))

        # pbc scroll bar
        self.pbc = QComboBox(); pbcBox.addWidget(self.pbc)
        self.pbc.addItems(['00', '11', '22'])
        self.pbc.currentIndexChanged[str].connect(self.selectPbc)
        self.string_widgets['pbc'] = self.pbc

        # Input part 2 - HBox. Second top region of the control part.
        inputLayout2 = QHBoxLayout()
        controlLayout.addLayout(inputLayout2)

        # Adding path text box to the layout
        pathBox = QVBoxLayout(); inputLayout2.addLayout(pathBox)
        pathBox.addWidget(QLabel("Input files path"))
        # Creating input files path text box 
        # In this case self is needed in order to reference later to the written text inside.
        self.path_inputFiles = QLineEdit(); pathBox.addWidget(self.path_inputFiles)
        self.path_inputFiles.setText('./')

        # Adding geometry_file text box to the layout
        geometryBox = QVBoxLayout(); inputLayout2.addLayout(geometryBox)
        geometryBox.addWidget(QLabel("Geometry file name"))
        # Creating geometry file path text box
        self.geometryFile = QLineEdit(); geometryBox.addWidget(self.geometryFile)
        self.geometryFile.setText('input.xyz')

        # Adding CP2K/GPAW name text box to the layout
        nameBox = QVBoxLayout(); inputLayout2.addLayout(nameBox)
        nameBox.addWidget(QLabel("CP2K/GPAW name"))
        # Creating name text box
        self.name = QLineEdit(); nameBox.addWidget(self.name)
        self.name.setText('none')

        # Button importing packages and everything
        importButton = QPushButton("Import"); inputLayout2.addWidget(importButton)
        importButton.clicked.connect(self.imported)

    #### Seperating line
        line = QFrame(); controlLayout.addWidget(line); line.setFrameShape(QFrame.HLine); line.setFrameShadow(QFrame.Sunken)

############################## GRID AND RUNNING OPTIONS ####################################

    ############# grl1 - grid running layout 1 - Format, K, Q ############################

        grl1 = QHBoxLayout(); controlLayout.addLayout(grl1); 

        # Adding format to layout
        grl1.addWidget(QLabel("         Format"))
        # Creating format scroll bar
        self.formatBar = QComboBox(); grl1.addWidget(self.formatBar);
        self.formatBar.addItems(['xsf', 'npy'])
        self.formatBar.currentIndexChanged[str].connect(self.selectFormat)
        self.string_widgets['data_format'] = self.formatBar

        # Adding K number range
        grl1.addWidget(QLabel('             K: '))
        self.k = QDoubleSpinBox(); grl1.addWidget(self.k);
        self.k.setRange(0.0, 2.0); self.k.setSingleStep(0.05)
        self.k.valueChanged.connect(self.selectK)
        self.k.setValue(0.24)
        self.num_widgets['kValue'] = self.k

        # Adding Q number range
        grl1.addWidget(QLabel('             Q: '))
        self.q = QDoubleSpinBox(); grl1.addWidget(self.q);
        self.q.setRange(0.0, 2.0); self.q.setSingleStep(0.05)
        self.q.valueChanged.connect(self.selectQ)
        self.num_widgets['qValue'] = self.q
        
    ############# grl2 - grid running layout 2 - Xmin, Xmax, dX #######################

        grl2 = QHBoxLayout(); controlLayout.addLayout(grl2);

        # Adding Xmin number range
        grl2.addWidget(QLabel('         Xmin: '))
        self.xMin = QDoubleSpinBox(self); grl2.addWidget(self.xMin);
        self.xMin.setRange(-50.0, 20.0); self.xMin.setSingleStep(0.05)
        self.xMin.valueChanged.connect(self.selectX)

        # Adding Xmax number range
        grl2.addWidget(QLabel('         Xmax: '))
        self.xMax = QDoubleSpinBox(); grl2.addWidget(self.xMax);
        self.xMax.setRange(0.0, 50.0); self.xMax.setSingleStep(0.05)
        self.xMax.setValue(20.0)
        self.xMax.valueChanged.connect(self.selectX)

        # Adding dX number range
        grl2.addWidget(QLabel('             dX: '))
        self.dx = QDoubleSpinBox(); grl2.addWidget(self.dx);
        self.dx.setRange(0.0, 20.0); self.dx.setSingleStep(0.05)
        self.dx.setValue(0.1)
        self.dx.valueChanged.connect(self.selectX)

    ############# grl3 - grid running layout 3 - Ymin, Ymax, dY ########################

        grl3 = QHBoxLayout(); controlLayout.addLayout(grl3);

        # Adding Ymin number range
        grl3.addWidget(QLabel('         Ymin: '))
        self.yMin = QDoubleSpinBox(); grl3.addWidget(self.yMin);
        self.yMin.setRange(-50.0, 20.0); self.yMin.setSingleStep(0.05)
        self.yMin.valueChanged.connect(self.selectY)

        # Adding Ymax number range
        grl3.addWidget(QLabel('         Ymax: '))
        self.yMax = QDoubleSpinBox(); grl3.addWidget(self.yMax);
        self.yMax.setRange(0.0, 50.0); self.yMax.setSingleStep(0.05)
        self.yMax.setValue(20.0)
        self.yMax.valueChanged.connect(self.selectY)

        # Adding dY number range
        grl3.addWidget(QLabel('             dY: '))
        self.dy = QDoubleSpinBox(); grl3.addWidget(self.dy);
        self.dy.setRange(0.0, 20.0); self.dy.setSingleStep(0.05)
        self.dy.setValue(0.1)
        self.dy.valueChanged.connect(self.selectY)

    ############# grl4 - grid running layout 4 - Zmin, Zmax, dZ #########################

        grl4 = QHBoxLayout(); controlLayout.addLayout(grl4);

        # Adding Zmin number range
        grl4.addWidget(QLabel('         Zmin: '))
        self.zMin = QDoubleSpinBox(); grl4.addWidget(self.zMin);
        self.zMin.setRange(-50.0, 20.0); self.zMin.setSingleStep(0.05)
        self.zMin.setValue(5.0)
        self.zMin.valueChanged.connect(self.selectZ)

        # Adding Zmax number range
        grl4.addWidget(QLabel('         Zmax: '))
        self.zMax = QDoubleSpinBox(); grl4.addWidget(self.zMax);
        self.zMax.setRange(0.0, 50.0); self.zMax.setSingleStep(0.05)
        self.zMax.setValue(6.0)
        self.zMax.valueChanged.connect(self.selectZ)

        # Adding dZ number range
        grl4.addWidget(QLabel('             dZ: '))
        self.dz = QDoubleSpinBox(); grl4.addWidget(self.dz);
        self.dz.setRange(0.0, 20.0); self.dz.setSingleStep(0.05)
        self.dz.setValue(1.0)
        self.dz.valueChanged.connect(self.selectZ)

    ############# grl5 - grid running layout 5  - Scan Type, Eta, WF_decay ##################

        grl5 = QHBoxLayout(); controlLayout.addLayout(grl5); 

        # Adding scan type to layout
        grl5.addWidget(QLabel("         Scan type"))
        # Creating scan type scroll bar
        self.scan = QComboBox(); grl5.addWidget(self.scan);
        self.scan.addItems(['dIdV', 'v-scan', 'states'])
        self.scan.currentIndexChanged[str].connect(self.selectScanType)
        self.string_widgets['scan_type'] = self.scan
        
        # Adding Eta number range
        grl5.addWidget(QLabel('             Eta:            '))
        self.eta = QDoubleSpinBox(); grl5.addWidget(self.eta);
        self.eta.setRange(0.0, 20.0); self.eta.setSingleStep(0.05)
        self.eta.valueChanged.connect(self.selectEta)
        self.eta.setValue(0.1)
        self.num_widgets['etaValue'] = self.eta

        # Adding WF_decay number range (wfd)
        grl5.addWidget(QLabel('             WF_decay: '))
        self.wfd = QDoubleSpinBox(); grl5.addWidget(self.wfd);
        self.wfd.setRange(0.0, 2.0); self.wfd.setSingleStep(0.05)
        self.wfd.valueChanged.connect(self.selectWF_decay)
        self.num_widgets['wf_decay'] = self.wfd

    ############# grl6 - grid running layout 6 - V(min), Vmax, dV ###########################

        grl6 = QHBoxLayout(); controlLayout.addLayout(grl6);

        # Adding V(min) number range
        grl6.addWidget(QLabel('         V(min): '))
        self.vMin = QDoubleSpinBox(); grl6.addWidget(self.vMin);
        self.vMin.setRange(-2.0, 2.0); self.vMin.setSingleStep(0.05)
        self.vMin.setValue(-2.0)
        self.num_widgets['V'] = self.vMin

        # Adding Vmax number range
        grl6.addWidget(QLabel('         Vmax: '))
        self.vMax = QDoubleSpinBox(); grl6.addWidget(self.vMax);
        self.vMax.setRange(0.0, 20.0); self.vMax.setSingleStep(0.05)
        self.vMax.setValue(2.0)
        self.num_widgets['Vmax'] = self.vMax

        # Adding dV number range
        grl6.addWidget(QLabel('             dV: '))
        self.dv = QDoubleSpinBox(); grl6.addWidget(self.dv);
        self.dv.setRange(0.0, 20.0); self.dv.setSingleStep(0.05)
        self.dv.setValue(0.1)
        self.num_widgets['dV'] = self.dv

        self.vMin.valueChanged.connect(self.selectV)
        self.dv.valueChanged.connect(self.selectV)
        self.vMax.valueChanged.connect(self.selectV)
        
    ############# grl7 - grid running layout 7 - Tip orb. s, Tip orb. pxy, OMP_NUM_THREADS ###########

        grl7 = QHBoxLayout(); controlLayout.addLayout(grl7);

        # Adding Tip orb s number range
        grl7.addWidget(QLabel('         Tip orb. s:  '))
        self.orbS = QDoubleSpinBox(); grl7.addWidget(self.orbS);
        self.orbS.setRange(-2.0, 2.0); self.orbS.setSingleStep(0.05)
        self.orbS.valueChanged.connect(self.selectTipOrbS)
        self.orbS.setValue(0.13)
        self.num_widgets['tipOrbS'] = self.orbS

        # Adding Vmax number range
        grl7.addWidget(QLabel('                 Tip orb. pxy:   '))
        self.orbPxy = QDoubleSpinBox(); grl7.addWidget(self.orbPxy);
        self.orbPxy.setRange(-2.0, 20.0); self.orbPxy.setSingleStep(0.05)
        self.orbPxy.valueChanged.connect(self.selectTipOrbPxy)
        self.orbPxy.setValue(0.87)
        self.num_widgets['tipOrbPxy'] = self.orbPxy

        # Adding dV number range
        grl7.addWidget(QLabel('             OMP_NUM_THREADS: '))
        self.ont = QSpinBox(); grl7.addWidget(self.ont);
        self.ont.setRange(0, 20); self.ont.setSingleStep(1)
        self.ont.valueChanged.connect(self.selectONT)
        self.ont.setValue(1)
        self.num_widgets['OMP_NUM_THREADS'] = self.ont

    ############# grl8 - grid running layout 8 - Tip Orb. , Tip type, Run button ###############

        grl8 = QHBoxLayout(); controlLayout.addLayout(grl8)

        l1 = QHBoxLayout(); grl8.addLayout(l1)

        # Tip type
        label2 = QLabel("               Tip type:")
        label2.setAlignment(Qt.AlignRight)
        l1.addWidget(label2)

        # Creating tip type scroll bar
        self.tip_type = QComboBox();  l1.addWidget(self.tip_type);
        self.tip_type.addItems(['fixed', 'relaxed'])
        self.tip_type.currentIndexChanged[str].connect(self.selectTipType)
        l1.setAlignment(Qt.AlignLeft)
        self.string_widgets['tip_type'] = self.tip_type

        # Buttons
        l2 = QHBoxLayout(); grl8.addLayout(l2)
        run = QPushButton("Run"); l2.addWidget(run)
        run.clicked.connect(self.running)

    #### Seperating line 
        line = QFrame(); controlLayout.addWidget(line); line.setFrameShape(QFrame.HLine); line.setFrameShadow(QFrame.Sunken)

############################# VISUALIZING OPTIONS PART ####################################

    ####### vo1 - visualizing options part 1, row 1 - Map type, Voltage, Height, Show button ############

        vo1 = QHBoxLayout(); controlLayout.addLayout(vo1)

        # Map type scroll bar
        vo1.addWidget(QLabel(       "Map type"))
        mapType = QComboBox(); vo1.addWidget(mapType)
        mapType.addItems(['dIdV', 'STM'])
        mapType.currentIndexChanged[str].connect(self.selectMapType)


        # Adding Voltage Index number range
        vo1.addWidget(QLabel('             Voltage Index: '))
        self.vIndx = QSpinBox(); vo1.addWidget(self.vIndx);
        self.vIndx.setRange(0, 10000); self.vIndx.setSingleStep(1)
        self.vIndx.valueChanged.connect(self.selectVindx)

        # Adding Height Index number range
        vo1.addWidget(QLabel('             Height Index: '))
        self.hIndx = QSpinBox(); vo1.addWidget(self.hIndx);
        self.hIndx.setRange(0, 10000); self.hIndx.setSingleStep(1)
        self.hIndx.valueChanged.connect(self.selectHindx)
        

        show = QPushButton("Show"); vo1.addWidget(show)
        show.clicked.connect(self.showImage)

    ####### vo2 - visualizing options part 2, row 2 - save image, save data, save options, load options (buttons) #################

        vo2 = QHBoxLayout(); controlLayout.addLayout(vo2)

        # Buttons
        saveImg = QPushButton("Save image"); vo2.addWidget(saveImg)
        saveImg.clicked.connect(self.saveFigure)

        saveData = QPushButton("Save data"); vo2.addWidget(saveData)
        saveData.clicked.connect(self.saveData)

        saveOptions = QPushButton("Save options"); vo2.addWidget(saveOptions)
        saveOptions.clicked.connect(self.saveOptions)

        loadOptions = QPushButton("Load options"); vo2.addWidget(loadOptions)
        loadOptions.clicked.connect(self.load)

        # Adding layout to the main widget 
        self.mainWidget.setFocus()
        self.setCentralWidget(self.mainWidget)

####################### Saving methods (figure, data, options) ####################################

    def saveFigure(self):
        if self.runClicked:
            self.myCanvas.saveFigure()
        else:
            print('Error, no plot exists. Run first.')
    
    def saveData(self):
        if self.runClicked:
            self.myCanvas.saveData()
        else:
            print('Error. Run first.')
    
    def saveOptions(self):
        fileName, fext = QFileDialog.getSaveFileName(self, "QFileDialog.getSaveFileName()", "", "ASCII File (*.txt)")
        if fileName:
            data = self.myDict
            paths = self.paths
            f = open(fileName, 'w')
            f.write("# comment can be written by hash and space afterwards \n")
            for d in data:
                f.write(d+': ')
                f.write(str(data[d]))
                f.write('\n')
            if paths:
                f.write('paths: ')
                f.write(str(paths))
            f.close()

############################# Select methods (selecting values for all vars)  ####################################

    def selectCode(self, code):
        self.myDict['dft_code'] = code

    def selectOrbital(self, orbital):
        self.myDict['sample_orbs'] = orbital

    def selectSpin(self, spin):
        self.myDict['spin'] = spin 

    def selectPbc(self, pbc):
        self.myDict['pbc'] = pbc
    
    def selectFormat(self, myFormat):
        self.myDict['data_format'] = myFormat

    def selectK(self):
        self.myDict['kValue'] = self.k.value()
    
    def selectQ(self):
        self.myDict['qValue'] = self.q.value()

    def selectX(self):
        self.myDict['x'][0] = self.xMin.value()
        self.myDict['x'][1] = self.xMax.value()
        self.myDict['x'][2] = self.dx.value()
    
    def selectY(self):
        self.myDict['y'][0] = self.yMin.value()
        self.myDict['y'][1] = self.yMax.value()
        self.myDict['y'][2] = self.dy.value()

    def selectZ(self):
        self.myDict['z'][0] = self.zMin.value()
        self.myDict['z'][1] = self.zMax.value()
        self.myDict['z'][2] = self.dz.value()

    def selectScanType(self, scanType):
        self.myDict['scan_type'] = scanType
    
    def selectEta(self):
        self.myDict['etaValue'] = self.eta.value()

    def selectWF_decay(self):
        self.myDict['wf_decay'] = self.wfd.value()

    def selectV(self):
        self.myDict['V'] = self.vMin.value()
        self.myDict['Vmax'] = self.vMax.value()
        self.myDict['dV'] = self.dv.value()
    
    def selectTipOrbS(self):
        self.myDict['tipOrbS'] = self.orbS.value()
    
    def selectTipOrbPxy(self):
        self.myDict['tipOrbPxy'] = self.orbPxy.value()
    
    def selectONT(self):
        self.myDict['OMP_NUM_THREADS'] = self.ont.value()
    
    
    def selectTipType(self, tip_type):
        self.myDict['tip_type'] = tip_type

    def selectMapType(self, chosenMap):
        self.map = chosenMap

    def selectVindx(self):
        self.Vindx = self.vIndx.value()
    
    def selectHindx(self):
        self.Hindx = self.hIndx.value()

############################# Button methods (import, load, run, show) ####################################

    def imported(self):
        # check if path to input files is given, if not give error and do not let to run
        inputPath = self.path_inputFiles.text()
        geometry_file = self.geometryFile.text() 
        cp2kName = self.name.text()

        if len(inputPath) == 0 or len(cp2kName) == 0 or len(geometry_file) == 0:
            print ('Path/name/file not given')
            return 

        self.paths = {'inputPath': inputPath, 
                      'geometry_file': geometry_file,
                      'cp2kName': cp2kName,}
        
        self.importData = importData(self.myDict, self.paths)

        if self.importData:
            print("energies prepared, coeffecients read")

    def load(self):
        filePath, ext = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "", "ASCII File (*.txt)")
        if filePath:
            print (" >> OVERWRITING SETTINGS by "+ filePath) 
            fin = open(filePath, 'r')
            for line in fin:
                line = re.sub('\[|\]|\:|\,|\)|\(|\'|\}|\{', '', line)
                words = line.split()
                if len(words) == 0 : continue
                if words[0][0] == '#': continue
                if words[0] in ['x', 'y', 'z']:
                    sampleList = [float(words[1]), float(words[2]), float(words[3])]
                    self.myDict[words[0]] = sampleList
                    self.paramList.remove(words[0])
                    continue
                if words[0] == 'paths':
                    self.paths = {}

                    # inputPath
                    self.paths[words[1]] = words[2]
                    self.path_inputFiles.setText(words[2])
                    
                    # geometry_file
                    self.paths[words[3]] = words[4]
                    self.geometryFile.setText(words[4])

                    # cp2kName
                    self.paths[words[5]] = words[6]
                    self.name.setText(words[6])

                    self.paramList.remove('paths')
                    continue

                if words[0] == 'pbc':
                    self.myDict[words[0]] = words[1]
                    self.paramList.remove(words[0])
                    continue
                # All the other cases. Assigning values to the keys of self.myDict and trying to convert them to float. If they are not a number conv2float will just return string
                self.myDict[words[0]] = conv2float(words[1])
                self.paramList.remove(words[0])

            self.updateWidgets()
            
            if len(self.paramList) != 0:
                print('Following parameters where not changed and remained default:')
                for p in self.paramList:
                    if p == 'paths':
                        print (p + ' = ' + str(self.paths))
                    else:
                        print(p + ' = ' + self.myDict[p])

            self.paramList = []            
            self.paramList = list(self.myDict.keys())
            self.paramList.append('paths')
            fin.close()

    def updateWidgets(self):

        self.xMin.blockSignals(True); self.xMin.setValue(self.myDict['x'][0]); self.xMin.blockSignals(False)
        self.xMax.blockSignals(True); self.xMax.setValue(self.myDict['x'][1]); self.xMax.blockSignals(False)
        self.dx.blockSignals(True); self.dx.setValue(self.myDict['x'][2]); self.dx.blockSignals(False)
        
        self.yMin.blockSignals(True); self.yMin.setValue(self.myDict['y'][0]); self.yMin.blockSignals(False)
        self.yMax.blockSignals(True); self.yMax.setValue(self.myDict['y'][1]); self.yMax.blockSignals(False)
        self.dy.blockSignals(True); self.dy.setValue(self.myDict['y'][2]); self.dy.blockSignals(False)

        self.zMin.blockSignals(True); self.zMin.setValue(self.myDict['z'][0]); self.zMin.blockSignals(False)
        self.zMax.blockSignals(True); self.zMax.setValue(self.myDict['z'][1]); self.zMax.blockSignals(False)
        self.dz.blockSignals(True); self.dz.setValue(self.myDict['z'][2]); self.dz.blockSignals(False)

        for ws, wn in zip(self.string_widgets, self.num_widgets):
            self.string_widgets[ws].blockSignals(True)
            self.string_widgets[ws].setCurrentIndex(self.string_widgets[ws].findText(self.myDict[ws]))
            self.string_widgets[ws].blockSignals(False)
            self.num_widgets[wn].blockSignals(True)
            self.num_widgets[wn].setValue(float(self.myDict[wn]))
            self.num_widgets[wn].blockSignals(False)
            
    def running(self):
        # run PPSTM_simple with either grid data from 'relaxed' tip or from 'fixed'
        if self.paths and self.importData:
            self.plotData = newPPSTM_simple(self.myDict, self.paths, self.importData)
        else:
            print ('Error. Not all the parameters were given (paths - click import button).')
            return 
        if self.plotData:
            # Update plot
            self.runClicked = True
            if self.plotImage(self.plotData, self.Hindx, self.Vindx, self.map):
                canvas = MplCanvas(parentApp=self, parentWidget=self.mainWidget, width=20, height=8, dpi=150,)
                self.layout.replaceWidget(self.myCanvas, canvas)
                self.myCanvas = canvas
    
    def showImage(self):
        if self.runClicked:
            if self.plotImage(self.plotData, self.Hindx, self.Vindx, self.map):
                canvas = MplCanvas(parentApp=self, parentWidget=self.mainWidget, width=20, height=8, dpi=150,)
                self.layout.replaceWidget(self.myCanvas, canvas)
                self.myCanvas = canvas
        else:
            print('Error. Run first.')
    
    def plotGeom(self, atomSize=0.1, edge=True, ec='k'):
        atoms, tmp1, tmp2 = Bu.loadAtoms(os.path.join(self.paths['inputPath'], 'input_plot.xyz')); del tmp1, tmp2;
        self.plotData['geom_plot'] = atoms
        plt.fig = plt.gcf()
        es = atoms[0]; xs = atoms[1]; ys = atoms[2]
        for i in range(len(xs)):
            fc = '#%02x%02x%02x' % elements.ELEMENT_DICT[es[i]][7] #; print "DEBUG: fc", fc ; ##fc = '#FFFFFF' ##
            if not edge:
                ec=fc
            circle=plt.Circle( ( xs[i], ys[i] ), atomSize, fc=fc, ec=ec  )
            plt.fig.gca().axes.add_artist(circle)
    
    def plotImage(self, plotData, height, voltage, mapType):

        if mapType == 'dIdV':
            if self.myDict['scan_type'] == 'states':
                print('Error. Change map type to STM')
                return None
            if height >= plotData['NoH_didv']:
                print ('Error. Height out of index')
                return None
            elif voltage >= plotData['NoV_didv']:
                print ('Error. Voltage out of index')
                return None
        elif mapType == 'STM':
            if self.myDict['scan_type'] == 'dIdV':
                print('Error. Change map type to dIdV')
                return None
            if height >= plotData['NoH_STM']:
                print ('Error. Height out of index')
                return None
            elif voltage >= plotData['NoV_STM']:
                print ('Error. Voltage out of index')
                return None
        print("We go to plotting ")
        # print "DEBUG: long name:::", namez[vv],';height:%03d;tip:'  %k,tip_type,';',tip_orb
        #name_plot=namez[voltage]+';height:'+str(height)+';tip:'+tip_type+';'+tip_orb
        plt.close()
        name_plot=plotData['namez'][voltage]+';height:'+str(height)+';tip:'+plotData['tip_type']
        if mapType == 'dIdV':
        # ploting part here:
            plt.figure( figsize=(0.5 * plotData['lvec'][1,0] , 0.5 * plotData['lvec'][2,1] ) )
            plt.imshow(plotData['didv'][voltage,height,:,:], origin='image', extent=plotData['extent'] , cmap='gray')
            self.plotGeom()
            plt.xlabel(r' Tip_x $\AA$')
            plt.ylabel(r' Tip_y $\AA$')
            plt.title("dIdV:"+name_plot)
        else:
        # ploting part here:
            plt.figure( figsize=(0.5 * plotData['lvec'][1,0] , 0.5 * plotData['lvec'][2,1] ) )
            plt.imshow(plotData['current'][voltage,height,:,:], origin='image', extent=plotData['extent'] , cmap='gray')
            self.plotGeom()
            plt.xlabel(r' Tip_x $\AA$')
            plt.ylabel(r' Tip_y $\AA$')
            plt.title("STM:"+name_plot)
    
        return True

if __name__ == "__main__":

    app = QApplication(sys.argv)

    window = Window()
    window.show()
    sys.exit(app.exec_())