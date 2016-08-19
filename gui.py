import os
import Tkinter
import Pmw

import chimera
from chimera.baseDialog import ModelessDialog
from chimera import openModels, Molecule, Element, Coord, selection
from Midas import wait
from chimera.colorTable import getColorByName
from chimera.colorTable import colors as colorDict
from chimera import runCommand as rc

import nmwiz
import numpy as np
import prody as pr
from StringIO import StringIO
import matplotlib.pyplot as plt

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class NMWizDialog(ModelessDialog):

	title = 'NMWiz - Main 1.2'
	name = 'Normal Mode Wizard - 1.2'

	buttons = ( 'Website', 'Close')
	help = 'prody.csb.pitt.edu'
	title = 'NMWiz - Main 1.2'
	Website = 'http://prody.csb.pitt.edu'

	def fillInUI(self, parent):

		self.file_name = Tkinter.StringVar(parent)
		LoadNMDbutton = Tkinter.Button(parent, text="Load NMD File", command= self.loadNMDFile)
		LoadNMDbutton.grid(column=0,row=0)
		FromMolbutton = Tkinter.Button(parent, text="From Molecule", command=self.FromMolecule)
		FromMolbutton.grid(column=0,row=1)
		FromMolbutton = Tkinter.Button(parent, text="ProDy Interface", command=self.prodyInterface)
		FromMolbutton.grid(column=0,row=2)
		FromMolbutton = Tkinter.Button(parent, text="Structure Comparison")
		FromMolbutton.grid(column=0,row=3)


	def loadNMDFile(self):

		self.parseNMDfile()
		if self.file_name.get() != "":
			self._buildNMDWindow()

			self._buildMolecule()
			self._showMolecule()
			self._buildArrows()
			self._drawArrows()

			self.molecules = openModels.list()
			self.proteinMol = self.molecules[-2]
			self.arrowMol = self.molecules[-1]
			rc("ribspline cardinal spec @CA")

			self._colorArrows()		

	def parseNMDfile(self):
		import tkFileDialog
		self.file_name.set(tkFileDialog.askopenfilename())

		if self.file_name.get() != "":
			length = file_len(self.file_name.get())
			k = 0
			mode_in = False
			f = open(self.file_name.get(), 'rb')
			for j in range(length):
				line = f.readline()
				prelabel = line.split(' ')[0]
				if prelabel == 'nmwiz_load':
					address = line.split(' ')[1]
					k += 1
				elif prelabel == 'name':
					self.name = line.split(' ')[1][:-1]
					k += 1
				elif prelabel == 'atomnames':
					self.atomNames = line.split(' ')[1:]
					self.atomNumber = len(self.atomNames)
					k += 1
				elif prelabel == 'resnames':
					self.resNames = line.split(' ')[1:]
					k += 1
				elif prelabel == 'resids':
					self.resIds = np.asarray(map(int,line.split(' ')[1:]))
					k += 1
				elif prelabel == 'chainids':
					self.chainIds = line.split(' ')[1:]
					k += 1
				elif prelabel == 'segnames':
					self.segNames = line.split(' ')[1:]
					k += 1
				elif prelabel == 'bfactors':
					self.bfactors = np.asarray(map(float,line.split(' ')[1:]))
					k += 1
				elif prelabel == 'coordinates':
					self.coordinates = np.asarray(map(float,line.split(' ')[1:])).reshape((self.atomNumber,3))
					k += 1
				elif prelabel == 'mode':
					if mode_in == False:
						modeNumber = length - k
						self.total_mode = modeNumber
						self.scales = np.zeros(modeNumber)
						self.scales2 = np.zeros(modeNumber)
						self.modeNo = np.zeros(modeNumber,dtype=int)
						self.modes = np.zeros((modeNumber,self.atomNumber*3))
						self.rmsd = Tkinter.StringVar()
						self.rmsd.set("2.0")
						mode_in = True
						mode_index = 0
					self.scales[mode_index] = float(line.split(' ')[2])
					self.modes[mode_index,:] = map(float,line.split(' ')[3:])
					self.scales2[mode_index] = float(self.rmsd.get())*np.sqrt(self.atomNumber)/self.scales[mode_index]
					self.modeNo[mode_index] = int(line.split(' ')[1])
					mode_index += 1
			
	def prodyInterface(self):
		self.prodyWindow = Tkinter.Toplevel()
		self.prodyWindow.title("NMWiz - ProDy Interface")

		group = Tkinter.LabelFrame(self.prodyWindow, text = "Atom Selection", padx = 5)
		group.pack()

		self.activeMoleculeStr = Tkinter.StringVar(self.prodyWindow)
		self.activeMoleculeId = Tkinter.IntVar(self.prodyWindow)
		self.atomNumberModel = Tkinter.IntVar(self.prodyWindow)
		self.selectionStr = Tkinter.StringVar(self.prodyWindow)
		self.selectionAtomNumber = Tkinter.IntVar(self.prodyWindow)
		molecules = openModels.list()
		moleculeNames = [m.name for m in molecules]
		moleculeIds = [m.id for m in molecules]
		self.molOptions = []
		molOptionsId = []
		for i in range(len(moleculeNames)):
			if isinstance(molecules[i],Molecule):
				self.molOptions.append(moleculeNames[i])
				molOptionsId.append(moleculeIds[i])

		if len(self.molOptions)<1:
			self.molOptions.append("")

		if len(moleculeNames)>0:
			self.activeMoleculeStr.set(self.molOptions[0])
			self.activeMoleculeId.set(0)
		else:
			self.activeMoleculeStr.set("")
			self.activeMoleculeId.set(0)

		w0_0 = Tkinter.Label(group, text="Molecule:")
		w0_0.grid(row=0,column=0)
		self.w0_1 = Tkinter.OptionMenu(group,self.activeMoleculeStr,()) #apply(Tkinter.OptionMenu, (group,self.active_mode) + tuple(OPTIONS))
		self.w0_1.grid(row=0,column=1,columnspan=3)
		self.w0_1['menu'].delete(0,'end')
		self.molecules = openModels.list()
		moleculeNames = [m.name for m in self.molecules]
		moleculeIds = [m.id for m in self.molecules]
		self.molOptions = []
		molOptionsId = []
		for i in range(len(moleculeNames)):
			if isinstance(self.molecules[i],Molecule):
				self.molOptions.append(moleculeNames[i])
				molOptionsId.append(moleculeIds[i])

		if len(self.molOptions)<1:
			self.molOptions = [""]

		for m in self.molOptions:
			self.w0_1['menu'].add_command(label=m,command=lambda value=m: self.updateActiveMolecule(value))

		w0_2 = Tkinter.Button(group, text="Update", command=self.updateMolecules)
		w0_2.grid(row=0,column=4,columnspan=2)
		self.w1_0 = Tkinter.Label(group, text="Information: " + str(self.atomNumberModel.get()) + " atoms")
		self.w1_0.grid(row=1,column=0,columnspan=6)
		if len(self.molOptions)>0:
			self.updateActiveMolecule(self.molOptions[0])

		w2_0 = Tkinter.Label(group, text="Selection:")
		w2_0.grid(row=2,column=0)
		self.selectionStr.set("@CA")
		self.selectionAtomNumber.set(0)
		w2_1 = Tkinter.Entry(group,textvariable=self.selectionStr)
		w2_1.grid(row=2,column=1,columnspan=3)
		w2_2 = Tkinter.Button(group, text="Select", command=self.selectMolecules)
		w2_2.grid(row=2,column=4,columnspan=2)


		self.w3_0 = Tkinter.Label(group, text="Information: " + str(self.selectionAtomNumber.get()) + " atoms are selected")
		self.w3_0.grid(row=3,column=0,columnspan=6)

		group2 = Tkinter.LabelFrame(self.prodyWindow, text = "Prody Job Settings", padx = 5)
		group2.pack()

		self.activeJob = Tkinter.StringVar(self.prodyWindow)
		self.outputDir = Tkinter.StringVar(self.prodyWindow)
		w4_0 = Tkinter.Label(group2, text="Prody Job:")
		w4_0.grid(row=0,column=0,columnspan=2)
		OPTIONS = ['ANM calculation','GNM calculation']
		self.activeJob.set(OPTIONS[0])
		w4_1 = Tkinter.OptionMenu(group2,self.activeJob,*(OPTIONS),command=self._changeActiveJob)
		w4_1.grid(row=0,column=2,columnspan=3)
		w5_0 = Tkinter.Label(group2, text="Output Directory:")
		w5_0.grid(row=1,column=0,columnspan=2)
		self.outputDir.set("")
		w5_1 = Tkinter.Entry(group2, textvariable=self.outputDir)
		w5_1.grid(row=1,column=2,columnspan=3)
		w5_2 = Tkinter.Button(group2, text="Browse", command=self.setOutputDir)
		w5_2.grid(row=1,column=6)
		w6_0 = Tkinter.Label(group2, text="Output Filename:")
		w6_0.grid(row=2,column=0,columnspan=2)
		self.outputFileName = Tkinter.StringVar(self.prodyWindow)
		w6_1 = Tkinter.Entry(group2, textvariable=self.outputFileName)
		w6_1.grid(row=2,column=2,columnspan=5)

		self.writeHeat = Tkinter.IntVar(self.prodyWindow)
		self.writeHeat.set(0)
		self.removeCoord = Tkinter.IntVar(self.prodyWindow)
		self.removeCoord.set(0)
		w7_0 = Tkinter.Checkbutton(group2, text="write and load cross-correlations heatmap", variable=self.writeHeat, onvalue=1, offvalue=0)
		w7_0.grid(row=3,column=0,columnspan=6)
		w8_0 = Tkinter.Checkbutton(group2, text="remove coordinate file upon job completion", variable=self.removeCoord, onvalue=1, offvalue=0)
		w8_0.grid(row=4,column=0,columnspan=6)

		group3 = Tkinter.LabelFrame(self.prodyWindow, text = "ANM Settings", padx = 5)
		group3.pack()

		self.numberOfModes = Tkinter.StringVar(self.prodyWindow)
		self.numberOfModes.set("10")
		w9_0 = Tkinter.Label(group3, text="Number of Modes:")
		w9_0.grid(row=0,column=0,columnspan=2)
		w9_1 = Tkinter.Entry(group3, textvariable=self.numberOfModes, width=6)
		w9_1.grid(row=0,column=2)
		self.Cutoff = Tkinter.StringVar(self.prodyWindow)
		self.Cutoff.set("15")
		w10_0 = Tkinter.Label(group3, text="Cutoff Distance(A):")
		w10_0.grid(row=1,column=0,columnspan=2)
		w10_1 = Tkinter.Entry(group3, textvariable=self.Cutoff, width=6)
		w10_1.grid(row=1,column=2)
		self.Gamma = Tkinter.StringVar(self.prodyWindow)
		self.Gamma.set("1")
		w10_2 = Tkinter.Label(group3, text="Force Constant:")
		w10_2.grid(row=1,column=3,columnspan=2)
		w10_3 = Tkinter.Entry(group3, textvariable=self.Gamma, width=6)
		w10_3.grid(row=1,column=5)

		# TO DO : Add extend model option

		w11_0 = Tkinter.Button(self.prodyWindow,text="Help",command=self._showHelpPrody,width=12)
		w11_0.pack(side=Tkinter.LEFT)
		w11_1 = Tkinter.Button(self.prodyWindow,text="Submit Job",command=self._submitJob,width=12)
		w11_1.pack(side=Tkinter.LEFT)
		w11_2 = Tkinter.Button(self.prodyWindow,text="ProDy Website", command=self._goWebSite,width=12)
		w11_2.pack(side=Tkinter.LEFT)

	def _goWebSite(self):
		import webbrowser
		webbrowser.open_new("http://prody.csb.pitt.edu")

	def _submitJob(self):
		try:
			import prody
			prody_imported = True
		except:
			prody_imported = False

		if prody_imported is True:
			coords = np.zeros((self.selectionAtomNumber.get(),3))
			atoms = selection.currentAtoms()
			coordOrder = []
			i = 0 
			for at in atoms:
				coordOrder.append(at.residue.id.position)
				for j in range(3):
					if j == 0:
						coords[i,j]=at.coord().x
					elif j == 1:
						coords[i,j]=at.coord().y
					elif j == 2:
						coords[i,j]=at.coord().z
				i += 1
			sortedCoord = np.argsort(np.array(coordOrder))
			coords = coords[sortedCoord,:]
			self.atomNumber = coords.shape[0]

			anm = prody.ANM(self.activeMoleculeStr.get())
			anm.buildHessian(coords,cutoff=float(self.Cutoff.get()),gamma=float(self.Gamma.get()))
			anm.calcModes(n_modes=int(self.numberOfModes.get()))
			atomGroup = prody.AtomGroup()
			atomGroup.setNames(['CA']*coords.shape[0])
			atomGroup.setCoords(coords)
			prody.writeNMD(self.outputDir.get() + "/" + self.outputFileName.get(), anm, atomGroup)
			self.total_mode = int(self.numberOfModes.get())

			self.scales = np.zeros(self.total_mode)
			self.scales2 = np.zeros(self.total_mode)
			self.modes = np.zeros((self.total_mode,coords.shape[0]*3))
			self.rmsd = Tkinter.StringVar()
			self.rmsd.set("2.0")
			eigvals = anm.getVariances()**0.5
			eigvecs = anm.getArray()
			for i in range(self.total_mode):
				self.scales[i] = eigvals[i]
				self.modes[i,:] = np.transpose(eigvecs[:,i])
				self.scales2[i] = float(self.rmsd.get())*np.sqrt(coords.shape[0])/self.scales[i]
			
			self.coordinates = coords
			print coords[0]
			self._buildNMDWindow()

			self._buildArrows()
			self._drawArrows()

			self.molecules = openModels.list()
			self.proteinMol = self.molecules[0]
			self.arrowMol = self.molecules[1]
			rc("ribspline cardinal spec @CA")

			self._colorArrows()		

	def _buildNMDWindow(self):

		self.arrow_direction = 1
		self.NMD_window = Tkinter.Toplevel()
		self.NMD_window.title("NMWiz - " + self.name)

		group = Tkinter.LabelFrame(self.NMD_window, text = self.name, padx = 5)
		group.pack()

		self.active_mode = Tkinter.StringVar(self.NMD_window)
		OPTIONS = list(map(str,range(1,self.total_mode+1)))
		self.active_mode.set(OPTIONS[0])

		w0_0 = Tkinter.Label(group, text="Active Mode:")
		w0_0.grid(row=0,column=0)
		w0_1 = Tkinter.OptionMenu(group,self.active_mode,*tuple(OPTIONS),command=self.updateMode) #apply(Tkinter.OptionMenu, (group,self.active_mode) + tuple(OPTIONS))
		w0_1.grid(row=0,column=1,columnspan=3)
		w0_2 = Tkinter.Button(group, text="<=", command=self._activeModeReduce)
		w0_2.grid(row=0,column=4)
		w0_3 = Tkinter.Button(group, text="+/-", command=self._arrowReverse)
		w0_3.grid(row=0,column=5)
		w0_4 = Tkinter.Button(group, text="=>", command=self._activeModeIncrease)
		w0_4.grid(row=0,column=6)

		self.activeColor = Tkinter.StringVar(self.NMD_window)
		COLOR_OPTIONS = colorDict.keys()
		self.activeColor.set(COLOR_OPTIONS[0])
		self.activeColorVariable = getColorByName(self.activeColor.get())

		w0_5 = Tkinter.OptionMenu(group,self.activeColor,*tuple(COLOR_OPTIONS),command=self._updateArrowColor)
		w0_5.grid(row=0,column=7,columnspan=3)
		
		self.current_scale = Tkinter.StringVar()
		self.current_scale2 = Tkinter.StringVar()
		w1_0 = Tkinter.Label(group, text="Scale by:")
		w1_0.grid(row=1,column=0)
		self.current_scale.set("%.2f" % self.scales[OPTIONS.index(self.active_mode.get())])
		w1_1 = Tkinter.Entry(group, textvariable=self.current_scale, width=4)
		w1_1.grid(row=1,column=1,columnspan=3)
		w1_2 = Tkinter.Label(group, text="x", width=1)
		w1_2.grid(row=1,column=4)
		self.current_scale2.set("%.2f" % self.scales2[OPTIONS.index(self.active_mode.get())])
		w1_3 = Tkinter.Entry(group, textvariable=self.current_scale2, width=6)
		w1_3.grid(row=1,column=5)
		w1_4 = Tkinter.Button(group, text="+1", command=self._incrementCurrentScale2)
		w1_4.grid(row=1,column=6)
		w1_5 = Tkinter.Button(group, text="+5", command=self._increment5CurrentScale2)
		w1_5.grid(row=1,column=7)
		w1_6 = Tkinter.Button(group, text="-5", command=self._decrement5CurrentScale2)
		w1_6.grid(row=1,column=8)
		w1_7 = Tkinter.Button(group, text="-1", command=self._decrementCurrentScale2)
		w1_7.grid(row=1,column=9)

		w2_0 = Tkinter.Label(group, text="RMSD (A):")
		w2_0.grid(row=2,column=0)
		w2_1 = Tkinter.Entry(group, textvariable=self.rmsd, width=4)
		w2_1.grid(row=2,column=1,columnspan=3)
		w2_2 = Tkinter.Button(group, text="+0.1", command=self._incrementRMSD)
		w2_2.grid(row=2,column=4)
		w2_3 = Tkinter.Button(group, text="+0.5", command=self._increment5RMSD)
		w2_3.grid(row=2,column=5)
		w2_4 = Tkinter.Button(group, text="-0.5", command=self._decrement5RMSD)
		w2_4.grid(row=2,column=6)
		w2_5 = Tkinter.Button(group, text="-0.1", command=self._decrementRMSD)
		w2_5.grid(row=2,column=7)

		w3_0 = Tkinter.Label(group, text="Selection:")
		w3_0.grid(row=3,column=0)
		self.selection = Tkinter.StringVar()
		self.selection.set("")
		w3_1 = Tkinter.Entry(group, textvariable=self.selection)
		w3_1.grid(row=3,column=1,columnspan=5)
		w3_2 = Tkinter.Button(group, text="Redraw", command=self._redraw)
		w3_2.grid(row=3,column=6,columnspan=2)

		w4_0 = Tkinter.Button(group, text="Plot Mobility", command=self._plotMobility)
		w4_0.grid(row=4,column=0,columnspan=3)
		w4_1 = Tkinter.Button(group, text="Load Heatmap", command=self._plotHeatmap)
		w4_1.grid(row=4,column=3,columnspan=3)
		w4_2 = Tkinter.Button(group, text="Plot Data", command=self._plotData)
		w4_2.grid(row=4,column=6,columnspan=3)

		w5_0 = Tkinter.Button(group, text="Main", command=self._showMain)
		w5_0.grid(row=5,column=0,columnspan=2)
		w5_1 = Tkinter.Button(group, text="Save", command=self._save)
		w5_1.grid(row=5,column=2,columnspan=2)
		w5_2 = Tkinter.Button(group, text="Remove", command=self._remove)
		w5_2.grid(row=5,column=4,columnspan=2)
		w5_3 = Tkinter.Button(group, text="Help", command=self._showHelp)
		w5_3.grid(row=5,column=6,columnspan=2)

		group2 = Tkinter.LabelFrame(self.NMD_window, text = "Actions", padx = 5)
		group2.pack()

		self.w6_0 = Tkinter.Label(group2, text="Mode (" + self.active_mode.get() + ")")
		self.w6_0.grid(row=0,column=0,columnspan=2)
		w6_1 = Tkinter.Button(group2, text="Draw", command=self._drawMode)
		w6_1.grid(row=0,column=2,columnspan=2)
		w6_2 = Tkinter.Button(group2, text="Clean", command=self._cleanMode)
		w6_2.grid(row=0,column=4,columnspan=2)
		self.w6_3 = Tkinter.Button(group2, text="Hide", command=self._hideMode)
		self.w6_3.grid(row=0,column=6,columnspan=2)
		self.w6_4 = Tkinter.Button(group2, text="Options", command=self._optionsMode)
		self.w6_4.grid(row=0,column=8,columnspan=2)

		w7_0 = Tkinter.Label(group2, text="Animation:")
		w7_0.grid(row=1,column=0,columnspan=2)
		w7_1 = Tkinter.Button(group2, text="Make", command=self._buildMovie)
		w7_1.grid(row=1,column=2,columnspan=2)
		self.w7_2 = Tkinter.Button(group2, text="Play", command=self._playMovie)
		self.w7_2.grid(row=1,column=4,columnspan=2)
		self.w7_3 = Tkinter.Button(group2, text="Hide", command=self._hideMovie)
		self.w7_3.grid(row=1,column=6,columnspan=2)
		self.w7_4 = Tkinter.Button(group2, text="Options", command=self._optionsMovie)
		self.w7_4.grid(row=1,column=8,columnspan=2)

		w8_0 = Tkinter.Label(group2, text="Figures:")
		w8_0.grid(row=2,column=0,columnspan=2)
		w8_1 = Tkinter.Button(group2, text="Clear", command=self._clearSelection)
		w8_1.grid(row=2,column=2,columnspan=2)
		w8_2 = Tkinter.Button(group2, text="Close", command=self._clearFigure)
		w8_2.grid(row=2,column=4,columnspan=2)
		self.w8_3 = Tkinter.Button(group2, text="Hide", command=self._hideSelection)
		self.w8_3.grid(row=2,column=6,columnspan=2)
		self.w8_4 = Tkinter.Button(group2, text="Options", command=self._optionsSelection)
		self.w8_4.grid(row=2,column=8,columnspan=2)

		w9_0 = Tkinter.Label(group2, text="Molecule:") ## Add the molecule number and check from nmwiz.tcl
		w9_0.grid(row=3,column=0,columnspan=2)
		w9_1 = Tkinter.Button(group2, text="Update", command=self._updateView)
		w9_1.grid(row=3,column=2,columnspan=2)
		w9_2 = Tkinter.Button(group2, text="Focus", command=self._focusView)
		w9_2.grid(row=3,column=4,columnspan=2)
		self.w9_3 = Tkinter.Button(group2, text="Hide", command=self._hideMolecule)
		self.w9_3.grid(row=3,column=6,columnspan=2)
		self.w9_4 = Tkinter.Button(group2, text="Options", command=self._optionsMolecule)
		self.w9_4.grid(row=3,column=8,columnspan=2)

	def setOutputDir(self):
		import tkFileDialog
		self.outputDir.set(tkFileDialog.askdirectory())

	def _changeActiveJob(self,active_job):
		self.activeJob.set(active_job)

	def selectMolecules(self):
		rc("select " + self.selectionStr.get())
		self.selectionAtomNumber.set(len(selection.currentAtoms()))
		self.w3_0.config(text="Information: " + str(self.selectionAtomNumber.get()) + " atoms are selected")



	def updateActiveMolecule(self,value):
		self.activeMoleculeStr.set(value)
		self.activeMoleculeId = self.molOptions.index(self.activeMoleculeStr.get())
		if self.activeMoleculeStr.get() != "":
			self.activeMolecule = self.molecules[self.activeMoleculeId]
			self.atomNumberModel.set(self.activeMolecule.atomCoordinatesArray().shape[0])
		else:
			self.activeMolecule = Molecule()
			self.atomNumberModel.set(0)
		self.w1_0.config(text="Information: " + str(self.atomNumberModel.get()) + " atoms")

	def updateMolecules(self):
		self.activeMoleculeStr.set("")
		self.w0_1['menu'].delete(0,'end')
		self.molecules = openModels.list()
		moleculeNames = [m.name for m in self.molecules]
		moleculeIds = [m.id for m in self.molecules]
		self.molOptions = []
		molOptionsId = []
		for i in range(len(moleculeNames)):
			if isinstance(self.molecules[i],Molecule):
				self.molOptions.append(moleculeNames[i])
				molOptionsId.append(moleculeIds[i])

		if len(self.molOptions)<1:
			self.molOptions = [""]
		for m in self.molOptions:
			self.w0_1['menu'].add_command(label=m,command=lambda value=m: self.updateActiveMolecule(value))

	def _buildMolecule(self):

		self.m = Molecule()
		for i in range(self.atomNumber):
			residues = []
			if i == 0:
				r = self.m.newResidue("ALA"," ",i," ")
				atomCA = self.m.newAtom("CA",Element("C"))
				atomCA.setCoord(Coord(self.coordinates[i,0],self.coordinates[i,1],self.coordinates[i,2]))
				r.addAtom(atomCA)
				atomCA_old = atomCA
				r.ribbonDisplay = True
				residues.append(r)
			else:
				r = self.m.newResidue("ALA"," ",i," ")
				atomCA = self.m.newAtom("CA",Element("C"))
				atomCA.setCoord(Coord(self.coordinates[i,0],self.coordinates[i,1],self.coordinates[i,2]))
				r.addAtom(atomCA)
				self.m.newBond(atomCA_old,atomCA)
				atomCA_old = atomCA
				r.ribbonDisplay = True
				residues.append(r)

	class NMWizTraj:
		def __len__(self):
			return len(self.molecule.coordSets)

		def __getitem__(self, i):
			return self.molecule.coordSets[i]

	def _clearSelection(self):
		rc("~select")

	def _clearFigure(self):
		plt.close()

	def _updateView():
		return

	def _focusView():
		rc("focus")

	def _buildMovie(self):
		self._createMovie()
		self._showMovie()

	def _createMovie(self):
		try:
			self.endFrame
		except:
			self.endFrame = 99
		try: 
			self.startFrame
		except:
			self.startFrame = 1
		self.ProteinMovie = Molecule()
		for i in range(self.atomNumber):
			residues = []
			if i == 0:
				r = self.ProteinMovie.newResidue("ALA"," ",i," ")
				atomCA = self.ProteinMovie.newAtom("CA",Element("C"))
				atomCA.setCoord(Coord(self.coordinates[i,0],self.coordinates[i,1],self.coordinates[i,2]))
				r.addAtom(atomCA)
				atomCA_old = atomCA
				r.ribbonDisplay = True
				residues.append(r)
			else:
				r = self.ProteinMovie.newResidue("ALA"," ",i," ")
				atomCA = self.ProteinMovie.newAtom("CA",Element("C"))
				atomCA.setCoord(Coord(self.coordinates[i,0],self.coordinates[i,1],self.coordinates[i,2]))
				r.addAtom(atomCA)
				self.ProteinMovie.newBond(atomCA_old,atomCA)
				atomCA_old = atomCA
				r.ribbonDisplay = True
				residues.append(r)

		crdSet1 = self.proteinMol.activeCoordSet
		natoms = self.atomNumber
		v = self.modes[int(self.active_mode.get())-1,:].reshape((self.atomNumber,3)) * float(self.current_scale2.get()) * float(self.current_scale.get())
		for f in range(self.startFrame, self.endFrame+1):
			setID = crdSet1.id + f 
			crdSet = self.ProteinMovie.newCoordSet(setID)
			m = Molecule()
			if f<25:
				for i in range(self.atomNumber):
					residues = []
					r = m.newResidue("ALA"," ",i," ")
					atomCA = m.newAtom("CA",Element("C"))
					atomCA.setCoord(Coord(self.coordinates[i,0]+v[i,0]*(f+1)*0.04,self.coordinates[i,1]+v[i,1]*(f+1)*0.04,self.coordinates[i,2]+v[i,2]*(f+1)*0.04))
					r.addAtom(atomCA)
					if i != 0:
						m.newBond(atomCA_old,atomCA)
					atomCA_old = atomCA
					r.ribbonDisplay = True
					residues.append(r)
			elif f<50:
				for i in range(self.atomNumber):
					residues = []
					r = m.newResidue("ALA"," ",i," ")
					atomCA = m.newAtom("CA",Element("C"))
					atomCA.setCoord(Coord(self.coordinates[i,0]+v[i,0]*(50-f-1)*0.04,self.coordinates[i,1]+v[i,1]*(50-f-1)*0.04,self.coordinates[i,2]+v[i,2]*(50-f-1)*0.04))
					r.addAtom(atomCA)
					if i != 0:
						m.newBond(atomCA_old,atomCA)
					atomCA_old = atomCA
					r.ribbonDisplay = True
					residues.append(r)

			elif f<75:
				for i in range(self.atomNumber):
					residues = []
					r = m.newResidue("ALA"," ",i," ")
					atomCA = m.newAtom("CA",Element("C"))
					atomCA.setCoord(Coord(self.coordinates[i,0]-v[i,0]*(f+1-50)*0.04,self.coordinates[i,1]-v[i,1]*(f+1-50)*0.04,self.coordinates[i,2]-v[i,2]*(f+1-50)*0.04))
					r.addAtom(atomCA)
					if i != 0:
						m.newBond(atomCA_old,atomCA)
					atomCA_old = atomCA
					r.ribbonDisplay = True
					residues.append(r)
			else:
				for i in range(self.atomNumber):
					residues = []
					r = m.newResidue("ALA"," ",i," ")
					atomCA = m.newAtom("CA",Element("C"))
					atomCA.setCoord(Coord(self.coordinates[i,0]-v[i,0]*(100-f-1)*0.04,self.coordinates[i,1]-v[i,1]*(100-f-1)*0.04,self.coordinates[i,2]-v[i,2]*(100-f-1)*0.04))
					r.addAtom(atomCA)
					if i != 0:
						m.newBond(atomCA_old,atomCA)
					atomCA_old = atomCA
					r.ribbonDisplay = True
					residues.append(r)

		 	for a1, a2 in zip(self.ProteinMovie.atoms, m.atoms):
		 		a1.setCoord(a2.coord(), crdSet)
		 	m.destroy()
		self.ensemble = self.NMWizTraj()
		self.ensemble.name = "NMWiz trajectory from NMWiz" 
		self.ensemble.startFrame = 1
		self.ensemble.endFrame = self.endFrame
		self.ensemble.molecule = self.ProteinMovie
	
	def _showMovie(self):	
		from Movie.gui import MovieDialog
		self.mov = MovieDialog(self.ensemble)
		rc("wait 1")
		rc("ribspline cardinal spec @CA")
		self.w7_3.config(text="Hide", command=self._hideMovie)

	def _playMovie(self):
		self.mov.startMovieForwardCallback()
		self.w7_2.config(text="Pause", command=self._pauseMovie)

	def _pauseMovie(self):
		self.mov.stopMovieCallback()
		self.w7_2.config(text="Play", command=self._playMovie)

	def _hideMovie(self):
		self.mov.Close()
		openModels.close(self.ProteinMovie)
		self._createMovie()
		self.w7_3.config(text="Show", command=self._showMovie)

	def _optionsMovie(self):
		self.group4 = Tkinter.LabelFrame(self.NMD_window, text = "Animation Options", padx = 5)
		self.group4.pack()
		self.w11_0 = Tkinter.Label(self.group4, text="Under Construction.")
		self.w11_0.grid(row=0, column=0)
		self.w7_4.config(command= self._alterOptionsMovie)

	def _alterOptionsMovie(self):
		self.w11_0.pack_forget()
		self.group4.pack_forget()
		self.w7_4.config(command= self._optionsMovie)	

	def _hideSelection(self):
		rc("~ribbon sel")
		self.w8_3.config(text="Show", command=self._showSelection)

	def _showSelection(self):
		rc("ribbon sel")
		self.w8_3.config(text="Hide", command=self._hideSelection)

	def _hideMolecule(self):
		openModels.close(self.proteinMol)
		self._buildMolecule()
		self.w9_3.config(text="Show", command=self._showMolecule)

	def _showMolecule(self):
		openModels.add([self.m])
		self.proteinMol = openModels.list()[0]
		rc("ribspline cardinal spec @CA")
		self.w9_3.config(text="Hide", command=self._hideMolecule)

	def _optionsMolecule(self):
		self.group6 = Tkinter.LabelFrame(self.NMD_window, text = "Molecule Options", padx = 5)
		self.group6.pack()
		self.w13_0 = Tkinter.Label(self.group6, text="Under Construction.")
		self.w13_0.grid(row=0, column=0)
		self.w9_4.config(command= self._alterOptionsMolecule)

	def _alterOptionsMolecule(self):
		self.w13_0.pack_forget()
		self.group6.pack_forget()
		self.w9_4.config(command= self._optionsMolecule)	

	def _optionsSelection(self):
		self.group5 = Tkinter.LabelFrame(self.NMD_window, text = "Figure Options", padx = 5)
		self.group5.pack()
		self.w12_0 = Tkinter.Label(self.group5, text="Under Construction.")
		self.w12_0.grid(row=0, column=0)
		self.w8_4.config(command= self._alterOptionsSelection)

	def _alterOptionsSelection(self):
		self.w12_0.pack_forget()
		self.group5.pack_forget()
		self.w8_4.config(command= self._optionsSelection)	


	def _calcMSF(self):
		v = self.modes[int(self.active_mode.get())-1,:].reshape((self.atomNumber,3))
		self.beta_active = np.sum(v*v,axis=1)

	def _calcCrossCorr(self):
		v = self.modes[int(self.active_mode.get())-1,:].reshape((self.atomNumber,3))
		self.crossCorr = np.dot(v,v.T)
	
	def _buildArrows(self):

		begin = self.coordinates
		v = self.modes[int(self.active_mode.get())-1,:].reshape((self.atomNumber,3))
		end = begin + v * float(self.current_scale2.get()) * float(self.current_scale.get())  
		stringDATA = ""
		for i in range(self.atomNumber):
			if np.linalg.norm(begin[i,:]-end[i,:])>1e-1:
				stringDATA = stringDATA + ".arrow %.3g %.3g %.3g %.3g %.3g %.3g %.3g %3.g\n" % (begin[i,0],begin[i,1],begin[i,2],end[i,0],end[i,1],end[i,2],0.2,0.4)
		self.arrows = StringIO(stringDATA)

	def _removeArrows(self):
		openModels.remove(self.arrowMol)

	def _drawArrows(self):
		self.arrowModel = openModels.open(self.arrows,type="Bild", identifyAs='Arrows ' + self.active_mode.get())
		self.arrowMol = openModels.list()[-1]

	def _colorArrows(self):
		self.arrowMol.color = getColorByName(self.activeColor.get())

	def _updateArrows(self):
		self._removeArrows()
		self._buildArrows()

	def _drawMode(self):
		self._drawArrows()
		self._colorArrows()

	def _cleanMode(self):
		self._updateArrows()

	def _hideMode(self):
		self._removeArrows()
		self._buildArrows()
		self.w6_3.config(text="Show", command=self._showMode)

	def _showMode(self):
		self._drawMode()
		self.w6_3.config(text="Hide", command=self._hideMode)

	def _optionsMode(self):
		self.group3 = Tkinter.LabelFrame(self.NMD_window, text = "Mode Graphics Options", padx = 5)
		self.group3.pack()
		self.w10_0 = Tkinter.Label(self.group3, text="Under Construction.")
		self.w10_0.grid(row=0, column=0)
		self.w6_4.config(command= self._alterOptionsMode)

	def _alterOptionsMode(self):
		self.w10_0.pack_forget()
		self.group3.pack_forget()
		self.w6_4.config(command= self._optionsMode)		

	def _updateLabel(self):
		self.w6_0.config(text="Mode (" + self.active_mode.get() + ")")

	def _plotMobility(self):
		self._calcMSF()
		plt.plot(self.resIds,self.beta_active,label="Mode " + self.active_mode.get())
		plt.xlabel("Atom/Residue #")
		plt.ylabel("")
		plt.title(self.name + " square fluctuations")
		legend = plt.legend()
		plt.show()

	def _plotHeatmap(self):
		self.heatMapFileName = Tkinter.StringVar()
		import tkFileDialog
		self.heatMapFileName.set(tkFileDialog.askopenfilename())
		if self.heatMapFileName.get() != "":
			heatmap, meta = nmwiz.parseHeatmap(heatmap = self.heatMapFileName.get())
		plt.contourf(self.resIds,self.resIds,heatmap)
		plt.xlabel("Atom/Residue #")
		plt.ylabel("Atom/Residue #")
		plt.title(meta['title'])
		plt.colorbar()
		plt.show()

	def _plotData(self):
		self.dataFileName = Tkinter.StringVar()
		import tkFileDialog
		self.dataFileName.set(tkFileDialog.askopenfilename())
		if self.dataFileName.get() != "":
			data = np.loadtxt(self.dataFileName.get())
		else:
			return
		if len(data)==self.atomNumber:
			plt.plot(self.resIds,data)
			plt.xlabel("Atom/Residue #")
			plt.ylabel("Data")
			plt.title("External Data")
			plt.show()
		else:
			raise ValueError('Size of data and size of protein do not match!')
		

	def _showMain(self):
		nmwiz_dialog()
		show_nmwiz_dialog()
		#self._toplevel.deiconify()
		
	def _save(self):
		import tkFileDialog
		f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".nmd")
		if f is None:
			return
		f.write("nmwiz_load " + f.name + "\n")
		f.write("name " + self.name + "\n")
		f.write("atomnames ")
		if len(self.atomNames) == self.atomNumber:
			for i in range(self.atomNumber):
				f.write(self.atomNames[i] + " ")
		f.write("\n")
		f.write("resnames ")
		if len(self.resNames) == self.atomNumber:
			for i in range(self.atomNumber):
				f.write(self.resNames[i] + " ")
		f.write("\n")
		f.write("resids ")
		if len(self.resIds) == self.atomNumber:
			for i in range(self.atomNumber):
				f.write(str(self.resIds[i]) + " ")
		f.write("\n")
		f.write("chainids ")
		if len(self.chainIds) == self.atomNumber:
			for i in range(self.atomNumber):
				f.write(self.chainIds[i] + " ")
		f.write("\n")
		f.write("bfactors ")
		if len(self.bfactors) == self.atomNumber:
			for i in range(self.atomNumber):
				f.write(str(self.bfactors[i]) + " ")
		f.write("\n")
		f.write("coordinates ")
		if len(self.coordinates) == self.atomNumber:
			for i in range(self.atomNumber*3):
				f.write(str(self.coordinates.reshape(self.atomNumber*3)[i]) + " ")
		f.write("\n")
		for i in range(self.total_mode):
			f.write("mode " + str(i+1) + " " + str(self.scales[i]) + " ")
			for j in range(self.atomNumber*3):
				f.write(str(self.modes[i,j]) + " ")
			f.write("\n")
		f.write("\n")
		f.close()

	def _remove(self):
		self.NMD_window.destroy()
		plt.close()
		openModels.remove(self.proteinMol)
		openModels.remove(self.arrowMol)

	def _showHelp(self):
		self.HelpWindow = Tkinter.Toplevel()
		self.HelpWindow.title("NMWiz Help")
		self.HelpText = Tkinter.Text(self.HelpWindow)
		self.HelpText.pack()
		self.HelpText.insert(Tkinter.END, nmwiz.helpTextWizard())

	def _showHelpPrody(self):
		self.HelpWindow = Tkinter.Toplevel()
		self.HelpWindow.title("NMWiz Help")
		self.HelpText = Tkinter.Text(self.HelpWindow)
		self.HelpText.pack()
		self.HelpText.insert(Tkinter.END, nmwiz.helpTextPrody())


		# print self.active_mode
	def _redraw(self):
		a = 1
		print a

	def _incrementCurrentScale2(self):
		self.current_scale2.set("%.2f" % (float(self.current_scale2.get()) + 1.0))
		self.rmsd.set("%.2f" % abs(float(self.current_scale2.get())*float(self.current_scale.get())/np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _increment5CurrentScale2(self):
		self.current_scale2.set("%.2f" % (float(self.current_scale2.get()) + 5.0))
		self.rmsd.set("%.2f" % abs(float(self.current_scale2.get())*float(self.current_scale.get())/np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _decrement5CurrentScale2(self):
		self.current_scale2.set("%.2f" % (float(self.current_scale2.get()) - 5.0))
		self.rmsd.set("%.2f" % abs(float(self.current_scale2.get())*float(self.current_scale.get())/np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _decrementCurrentScale2(self):
		self.current_scale2.set("%.2f" % (float(self.current_scale2.get()) - 1.0))
		self.rmsd.set("%.2f" % abs(float(self.current_scale2.get())*float(self.current_scale.get())/np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _incrementRMSD(self):
		self.rmsd.set("%.2f" % (float(self.rmsd.get()) + 0.1))
		self.current_scale2.set("%.2f" % abs(float(self.rmsd.get())/float(self.current_scale.get())*np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _increment5RMSD(self):
		self.rmsd.set("%.2f" % (float(self.rmsd.get()) + 0.5))
		self.current_scale2.set("%.2f" % abs(float(self.rmsd.get())/float(self.current_scale.get())*np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _decrement5RMSD(self):
		self.rmsd.set("%.2f" % (float(self.rmsd.get()) - 0.5))
		self.current_scale2.set("%.2f" % abs(float(self.rmsd.get())/float(self.current_scale.get())*np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def _decrementRMSD(self):
		self.rmsd.set("%.2f" % (float(self.rmsd.get()) - 0.1))
		self.current_scale2.set("%.2f" % abs(float(self.rmsd.get())/float(self.current_scale.get())*np.sqrt(self.atomNumber)))
		self._updateArrowScale()

	def updateMode(self,active_mode):
		OPTIONS = list(map(str,range(1,self.total_mode+1)))
		COLOR_OPTIONS = colorDict.keys()
		self.current_scale.set("%.2f" % self.scales[OPTIONS.index(active_mode)])
		self.current_scale2.set("%.2f" % self.scales2[OPTIONS.index(active_mode)])
		self.activeColor.set(COLOR_OPTIONS[OPTIONS.index(active_mode)])
		self._updateArrows()
		self._drawArrows()
		self._colorArrows()
		self._updateLabel()

	def _updateArrowScale(self):
		self._updateArrows()
		self._drawArrows()

	def _updateArrowColor(self,arrow_color):
		self.arrowMol.color = getColorByName(arrow_color)
		self._colorArrows()

	def _activeModeReduce(self):
		c = (int(self.active_mode.get())-1)%self.total_mode
		if c == 0:
			self.active_mode.set(self.total_mode)
		else:
			self.active_mode.set(c)
		self.updateMode(self.active_mode.get())
		
	def _arrowReverse(self):
		self.arrow_direction *= -1
		self.current_scale2.set("%.2f" % (float(self.current_scale2.get())*-1))
		self._updateArrows()
		self._drawArrows()

	def _activeModeIncrease(self):
		self.active_mode.set((int(self.active_mode.get())+1)%self.total_mode)
		self.updateMode(self.active_mode.get())
		
	def FromMolecule(self):
		print self.writeHeat.get()

	def Show(self):

		nmwiz.showNMD(atomModeMap[atomMode.get()],
		bondModeMap[bondMode.get()])



def nmwiz_dialog(create = False):

	from chimera import dialogs
	return dialogs.find(NMWizDialog.name, create=create)
	
# -----------------------------------------------------------------------------
#
def show_nmwiz_dialog():

	from chimera import dialogs
	return dialogs.display(NMWizDialog.name)

# -----------------------------------------------------------------------------
#
from chimera import dialogs
dialogs.register(NMWizDialog.name, NMWizDialog, replace = True)
