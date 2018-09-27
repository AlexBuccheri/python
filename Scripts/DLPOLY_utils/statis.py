#!/usr/bin/env python3 
# -*- coding: UTF-8 -*-

__author__ = "Alin Marin Elena <alin@elena.space>"
__copyright__ = "Copyright© 2018 Alin M Elena"
__license__ = "GPL-3.0-only"
__version__ = "1.0"

import sys
from PyQt5.QtWidgets import QWidget,QApplication,QMainWindow,QLineEdit,\
        QPushButton,QLabel,QAction,QTableWidget,QTableWidgetItem,QVBoxLayout,\
        QHBoxLayout,QTabWidget,QGroupBox,QFormLayout,QComboBox,QSpinBox,\
        QDoubleSpinBox,QRadioButton,QSizePolicy,QCheckBox,QGridLayout

from PyQt5.QtGui import QIcon,QIntValidator,QDoubleValidator
from PyQt5.QtCore import pyqtSlot, Qt
if hasattr(Qt, 'AA_EnableHighDpiScaling'):
  QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    
if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
  QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure
import numpy as np
    

countries = ['A','B','C']

class App(QMainWindow):
  def __init__(self):
    super().__init__()
    self.title="DL_POLY statis file explorer"
    self.left=0
    self.top=0
    s = app.primaryScreen().size()
    self.width=s.width()
    self.height=s.height()
    self.InitUI()

  def InitUI(self):
    self.setWindowTitle(self.title)
    self.setGeometry(self.left,self.top,self.width,self.height)
    self.statusBar().showMessage("earth, milky way")
    
    self.main = MyTabs(self)
    self.setCentralWidget(self.main)
    #self.createActions()
    #self.main.layout.addWidget(self.act)

    menu = self.menuBar()
    fileMenu = menu.addMenu('File')
    helpMenu = menu.addMenu('Help')

#  file menu buttons
    exit = QAction(QIcon('exit.png'), 'Exit', self)
    exit.setShortcut('Ctrl+Q')
    exit.setStatusTip('Exit application')
    exit.triggered.connect(self.close)
    fileMenu.addAction(exit)

    exit = QAction(QIcon('save.png'), 'Flatten STATIS', self)
    exit.setShortcut('Ctrl+F')
    exit.setStatusTip('Flatten Statis file in timeseries')
    exit.triggered.connect(self.flatten)
    fileMenu.addAction(exit)
    
    self.show()

  def flatten(self):
    n,nd,data,names = readStatis(filename="STATIS")
    for i in range(nd-3):
      with open(names[i],'w') as f:
        for j in range(n):
            f.write("{0} {1}\n".format(data[j,1],data[j,i+3]))

class MyTabs(QWidget):
  def __init__(self,parent):
    super(QWidget, self).__init__(parent)

    self.layout = QVBoxLayout()

    self.tabs = QTabWidget() 
    self.tab_statis = QWidget()
    self.tab_rdf = QWidget()
    self.tab_advanced = QWidget()

    self.tabs.addTab(self.tab_statis,"STATIS")

    self.tab_statis.layout = QVBoxLayout(self)
    self.plot = FigureCanvas(Figure(figsize=(16,12)))
    self.tab_statis.layout.addWidget(NavigationToolbar(self.plot, self))
    self.tab_statis.layout.addWidget(self.plot)
    self.tab_statis.layout.addStretch(1) 
    self.n,self.nd,self.data,self.datumNames = readStatis()
    self.createTimeSeries(self.datumNames)
    self.tab_statis.layout.addWidget(self.ts)
    self.ax = self.plot.figure.add_subplot(111)
    self.tsb.currentIndexChanged.connect(self.load_data)
    self.tsb.currentIndexChanged.connect(self.update_limits)
    self.update_plot(self.tsb.currentIndex())
    self.createAnalysis()
    self.tab_statis.layout.addWidget(self.an)
    self.ana.currentIndexChanged.connect(self.do_analysis)
    self.bins.editingFinished.connect(self.do_histogram)
    self.pdf.stateChanged.connect(self.do_norm)
    self.lag.editingFinished.connect(self.do_autocorr)
    self.nsample.editingFinished.connect(self.do_run_ave)
    self.xmin.editingFinished.connect(self.do_run_ave)
    self.xmax.editingFinished.connect(self.do_run_ave)
    self.tab_statis.setLayout(self.tab_statis.layout)

    self.tabs.addTab(self.tab_advanced,"Advanced")
    self.tab_advanced.layout = QVBoxLayout(self)
    self.lout = QLabel("Output")
    self.tab_advanced.layout.addWidget(self.lout)
    self.createFormGen()
    self.tab_advanced.layout.addWidget(self.fgen)
    self.tab_advanced.setLayout(self.tab_advanced.layout)

    self.tabs.addTab(self.tab_rdf,"RDFDAT")
    self.tab_rdf.layout = QVBoxLayout(self)
    self.rdf_plot = FigureCanvas(Figure(figsize=(16,12)))
    self.tab_rdf.layout.addWidget(NavigationToolbar(self.rdf_plot, self))
    self.tab_rdf.layout.addWidget(self.rdf_plot)
    self.tab_rdf.layout.addStretch(1) 
    self.nrdf,self.nprdf,self.rdf_data,self.rdf_labels=readRDF(filename="RDFDAT")
    self.axr = self.rdf_plot.figure.add_subplot(111)
    if self.nrdf>0 :
      self.createRDFS(self.rdf_labels)
      self.tab_rdf.layout.addWidget(self.rg)
      self.tab_rdf.layout.addStretch(1)
      self.rdfs.currentIndexChanged.connect(self.load_rdf)
      self.update_rdf(self.rdfs.currentIndex())
    self.tab_rdf.setLayout(self.tab_rdf.layout)

    self.layout.addWidget(self.tabs)
    self.setLayout(self.layout)

  def createRDFS(self,rdfnames):
    layout=QHBoxLayout()
    self.rg=QGroupBox("Data")
    self.rdfs = QComboBox()
    self.rdfs.addItems(rdfnames)
    self.rdfs.setCurrentIndex(0)
    label=QLabel("Data set")
    layout.addWidget(label)
    layout.addWidget(self.rdfs)
    layout.addStretch(1)
    self.rg.setLayout(layout)

  def createTimeSeries(self,datumNames):
      self.ts = QGroupBox("Data")
      layout = QHBoxLayout()
      self.tsb = QComboBox()
      self.tsb.addItems(datumNames)
      self.tsb.setCurrentIndex(0)
      label = QLabel("Data Set")
      layout.addWidget(label)
      layout.addWidget(self.tsb)
      layout.addStretch(1)
      self.ts.setLayout(layout)

  def createAnalysis(self):
      self.an = QGroupBox("Toobox")
      layout = QHBoxLayout()
      self.ana = QComboBox()
      self.ana.addItems(['Timeseries','Histogram','Running average','Autocorrelation','Fourier Transform'])
      self.ana.setCurrentIndex(0)
      label = QLabel("Analysis:")
      layout.addWidget(label)
      layout.addWidget(self.ana)
      self.do_options()
      layout.addWidget(self.tool_options)
      self.an.setLayout(layout)

  def do_analysis(self,c):
      self.toggle_ana_options(c)
      if c == 0:
          self.update_plot(self.tsb.currentIndex())
      elif c == 1:    
          self.do_histogram(self.tsb.currentIndex())
      elif c == 2:    
          self.do_run_ave(self.tsb.currentIndex())
      elif c == 3:    
          self.do_autocorr(self.tsb.currentIndex())
      elif c == 4:    
          self.do_ft(self.tsb.currentIndex())

  def load_data(self,c):
      self.do_analysis(self.ana.currentIndex())

  def load_rdf(self,c):
      
    self.update_rdf(self.rdfs.currentIndex())

  def update_rdf(self,c):
    self.axr.clear()
    self.axr.plot(self.rdf_data[c,:,0],self.rdf_data[c,:,1], 'r-')
    self.axr.set_title(self.rdf_labels[c])
    self.axr.set_xlabel('r [A]')
    self.axr.set_ylabel('check units')
    self.axr.set_xlim(left=0.0,right=max(self.rdf_data[c,:,0]))
    self.rdf_plot.draw()

    
  def update_plot(self,c):
    self.ax.clear()
    self.ax.plot(self.data[:,1],self.data[:,c+3], 'r-')
    self.ax.set_title(self.datumNames[c])
    self.ax.set_xlabel('Time [ps]')
    self.ax.set_ylabel('check units')
    self.ax.set_xlim(left=min(self.data[:,1]),right=max(self.data[:,1]))
    self.plot.draw()

  def do_run_ave(self,c=-1):
    c= self.tsb.currentIndex()
    self.ax.clear()
    xmin=float(self.xmin.text())
    xmax=float(self.xmax.text())
    y = self.data[:,1]
    inx=np.intersect1d(np.where(y>=xmin),np.where(y<=xmax))[::int(self.nsample.text())]
    x =np.empty([len(inx)],dtype=float)
    y =np.empty([len(inx)],dtype=float)
    v =np.empty([len(inx)],dtype=float)
    s=0.0
    k=0
    for i in inx:
      k=k+1
      s = s + self.data[i,c+3]
      y[k-1] = self.data[i,c+3]
      v[k-1] = s/k
      x[k-1] = self.data[i,1]
    av = v[k-1] 
    ic=self.error.currentIndex()
    se=0.0
    if ic == 0:
      se = np.std(y)
    elif ic == 1:
      w=int(self.nwin.text())
      nb=len(y)//w
      avj=[np.average(y[i*w:(i+1)*w]) for i in range(nb)]
      var=np.sum((avj-av)*(avj-av)) 
      se = np.sqrt(var/nb/(nb-1))
    elif ic == 2:
      w=int(self.nwin.text())
      n=len(y)
      nb=n//w
      avj=[np.average(y[i*w:(i+1)*w]) for i in range(nb)]
      nav=[ (n*av-w*avj[i])/(n-w) for i in range(nb) ]
      var=np.sum((nav-av)*(nav-av))
      se = np.sqrt(var/nb*(nb-1))


    self.ave.setText("{0:.8E} +/- {1:.8E}".format(v[k-1],se))  
    self.ax.plot(self.data[:,1],self.data[:,c+3], 'r-',x,v,'b-')
    self.ax.set_title(self.datumNames[c])
    self.ax.set_xlabel('Time [ps]')
    self.ax.set_ylabel('check units')
    self.ax.set_xlim(left=min(self.data[:,1]),right=max(self.data[:,1]))
    self.plot.draw()


  def do_autocorr(self,l=-1):
    self.ax.clear()
    l = int(self.lag.text())
    c= self.tsb.currentIndex()
    av = np.average(self.data[:,c+3])
    var = np.var(self.data[:,c+3])
    y=self.data[:,c+3]-av
    n=len(y)
    x=np.arange(0,l,1)
    ci=np.correlate(y,y,mode='full')[-n:]
    ci=np.array([(y[:n-k]*y[-(n-k):]).sum() for k in range(l)])
    ci=ci/(var*(np.arange(n,n-l,-1)-x))
    
    self.ax.plot(x,ci, 'r-')
    self.ax.set_title(self.datumNames[c])
    self.ax.set_xlabel('Lag')
    self.ax.set_ylabel('')
    self.ax.set_xlim(left=0,right=l)
    self.plot.draw()

  def do_norm(self,c):
      self.do_histogram()
  
  def do_histogram(self,c=-1):
    if c == -1 :
        c= self.tsb.currentIndex()
    self.ax.clear()
    hist,bins = np.histogram(self.data[:,c+3],bins=int(self.bins.text()),density=self.pdf.isChecked())
    self.ax.set_title(self.datumNames[c])
    self.ax.set_xlabel('bins')
    self.ax.set_ylabel('pdf')
    self.ax.bar(bins[:-1],hist)
    self.plot.draw()

  def do_ft(self,c):
    self.ax.clear()
    ft = np.fft.rfft(self.data[:,c+3])
    d=self.data[1,1]-self.data[0,1]
    freq = np.fft.rfftfreq(self.n, d=d)
    self.ax.set_xlabel('Frequency [Hz]')
    self.ax.set_ylabel('')
    self.ax.set_title("FT -"+self.datumNames[c])
    self.ax.plot(freq,2*np.abs(ft)/self.n,'g-')
    self.plot.draw()

  def update_limits(self,c):
    self.xmax.setText('{0:10.4f}'.format(self.data[-1,1]))
    self.xmin.setText('{0:10.4f}'.format(self.data[0,1]))

  def do_options(self):    
    self.tool_options = QGroupBox('Options')
    lt = QHBoxLayout()

    self.hist_options = QGroupBox("Histogram")
    layout = QFormLayout()
    self.bins = QLineEdit()
    self.bins.setValidator(QIntValidator())
    self.bins.setText('42')
    layout.addRow(QLabel('Bins '),self.bins)
    self.pdf = QCheckBox()
    layout.addRow(QLabel('Normalised? '),self.pdf)
    self.hist_options.setLayout(layout)


    self.ave_options = QGroupBox("Average")
    layout = QFormLayout()
    self.nsample = QLineEdit()
    self.nsample.setValidator(QIntValidator())
    self.nsample.setText('1')
    layout.addRow(QLabel('Sample '),self.nsample)
    self.xmin = QLineEdit()
    self.xmin.setValidator(QDoubleValidator())
    self.xmin.setText('{0:10.4f}'.format(self.data[0][1]))
    layout.addRow(QLabel('xmin'),self.xmin)
    self.xmax = QLineEdit()
    self.xmax.setValidator(QDoubleValidator())
    self.xmax.setText('{0:10.4f}'.format(self.data[self.n-1][1]))
    layout.addRow(QLabel('xmax'),self.xmax)
    self.error = QComboBox()
    self.error.addItems(['std error','binning','jackknife'])
    self.error.setCurrentIndex(0)
    self.error.currentIndexChanged.connect(self.do_error)
    lh=QHBoxLayout()
    lh.addWidget(QLabel('Error'))
    lh.addWidget(self.error)
    lh.addWidget(QLabel('Window size: '))
    self.nwin = QLineEdit()
    self.nwin.setValidator(QIntValidator())
    self.nwin.setText('100')
    self.nwin.setEnabled(False)
    lh.addWidget(self.nwin)
    layout.addRow(lh)
    self.ave =QLabel()
    layout.addRow(QLabel('Average:  '),self.ave)
    self.ave_options.setLayout(layout)

    self.auto_options = QGroupBox("AutoCorrelation")
    layout = QFormLayout()
    self.lag = QLineEdit()
    self.lag.setValidator(QIntValidator())
    self.lag.setText('250')
    layout.addRow(QLabel('Lag '),self.lag)
    self.auto_options.setLayout(layout)


    self.toggle_ana_options(0)
    lt.addWidget(self.hist_options)
    lt.addWidget(self.ave_options)
    lt.addWidget(self.auto_options)
    self.tool_options.setLayout(lt)

  def do_error(self,c):
      self.nwin.setEnabled(False)
      if c>0:
          self.nwin.setEnabled(True)
      self.do_run_ave()

  def toggle_ana_options(self,c):
    self.auto_options.setEnabled(False)
    self.ave_options.setEnabled(False)
    self.hist_options.setEnabled(False)
    if c == 1:
        self.hist_options.setEnabled(True)
    if c == 2:
        self.ave_options.setEnabled(True)
    if c == 3:
        self.auto_options.setEnabled(True)


  @pyqtSlot()
  def on_click(self):
      b_col = self.sender()
      print("colour selected {0:s}".format(b_col.text()))

  def readForm(self):
    print(self.name.text())
    print(self.age.text())
    print(self.country.currentText())
    print(self.height.text())

  def createFormGen(self):

    self.fgen = QGroupBox("AA")
    layout = QFormLayout()
    self.name = QLineEdit()
    self.name.setText('John')
    self.name.setFixedWidth(200)
    layout.addRow(QLabel('Name'),self.name)
    self.country = QComboBox()
    self.country.addItems(countries)
    layout.addRow(QLabel('Country'),self.country)
    self.age = QLineEdit()
    self.age.setValidator(QIntValidator())
    self.age.setText('42')
    layout.addRow(QLabel('Age'),self.age)
    self.height = QLineEdit()
    self.height.setValidator(QDoubleValidator())
    self.height.setText('42.0')
    layout.addRow(QLabel('Height'),self.height)
    self.fgen.setLayout(layout)

def readRDF(filename="RDFDAT"):
    try:
      title, header, rdfall = open(filename).read().split('\n',2)
    except IOError:
        return 0,0,0,[]
    nrdf,npoints=map(int, header.split())
    b=2*(npoints+1)
    d=np.zeros((nrdf,npoints,2),dtype=float)
    labels=[]
    s=rdfall.split()
    for i in range(nrdf):
      x=s[b*i:b*(i+1)]
      y=np.array(x[2:],dtype=float)
      y.shape= npoints,2
      d[i,:,:]=y
      labels.append(x[0]+" ... "+x[1])
    return nrdf,npoints,d,labels

def readStatis(filename="STATIS"):
  h1, h2, s = open(filename).read().split('\n',2)
  d = np.array(s.split(), dtype=float)
  nd = int(d[2])
  n = d.size//(nd+3)
  d.shape = n, nd+3
  datumNames = [ "1-1 total extended system energy",
             "1-2 system temperature",
             "1-3 configurational energy",
             "1-4 short range potential energy",
             "1-5 electrostatic energy",
             "2-1 chemical bond energy",
             "2-2 valence angle and 3-body potential energy",
             "2-3 dihedral, inversion, and 4-body potential energy",
             "2-4 tethering energy",
             "2-5 enthalpy (total energy + PV)",
             "3-1 rotational temperature",
             "3-2 total virial",
             "3-3 short-range virial",
             "3-4 electrostatic virial",
             "3-5 bond virial",
             "4-1 valence angle and 3-body virial",
             "4-2 constraint bond virial",
             "4-3 tethering virial",
             "4-4 volume",
             "4-5 core-shell temperature",
             "5-1 core-shell potential energy",
             "5-2 core-shell virial",
             "5-3 MD cell angle α",
             "5-4 MD cell angle β",
             "5-5 MD cell angle γ",
             "6-1 PMF constraint virial",
             "6-2 pressure",
             "6-3 exdof"]

  for i in range(28,nd):
      datumNames.append("{0:d}-{1:d} col_{2:d}".format(i//5+1,i%5+1,i+1))
  return n,nd,d,datumNames

if __name__ == '__main__':
  app = QApplication(sys.argv)
  exe = App()
  sys.exit(app.exec_())
