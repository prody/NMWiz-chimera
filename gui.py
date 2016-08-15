import os
import Tkinter
import Pmw

import chimera
from chimera.baseDialog import ModelessDialog
from chimera import openModels, Molecule, Element, Coord
from Midas import wait
from chimera.colorTable import getColorByName
from chimera.colorTable import colors as colorDict

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
		FromMolbutton = Tkinter.Button(parent, text="ProDy Interface")
		FromMolbutton.grid(column=0,row=2)
		FromMolbutton = Tkinter.Button(parent, text="Structure Comparison")
		FromMolbutton.grid(column=0,row=3)


	def loadNMDFile(self):
		import tkFileDialog
		self.file_name.set(tkFileDialog.askopenfilename())

		if self.file_name.get() != "":
			length = file_len(self.file_name.get())
			f = open(self.file_name.get(), 'rb')
			line = f.readline()
			address = line 
			line = f.readline()
			self.name = line.split(' ')[1][:-1]
			line = f.readline()
			self.atomNames = line.split(' ')[1:]
			self.atomNumber = len(self.atomNames)
			line = f.readline()
			self.resNames = line.split(' ')[1:]
			line = f.readline()
			self.resIds = np.asarray(map(int,line.split(' ')[1:]))
			line = f.readline()
			self.chainIds = line.split(' ')[1:]
			line = f.readline()
			line = f.readline()
			self.bfactors = np.asarray(map(float,line.split(' ')[1:]))
			line = f.readline()
			self.coordinates = np.asarray(map(float,line.split(' ')[1:])).reshape((self.atomNumber,3))
			modeNumber = length - 9
			self.total_mode = modeNumber
			self.scales = np.zeros(modeNumber)
			self.scales2 = np.zeros(modeNumber)
			self.modes = np.zeros((modeNumber,self.atomNumber*3))
			self.rmsd = Tkinter.StringVar()
			self.rmsd.set("2.0")
			for i in range(modeNumber):
				line = f.readline()
				self.scales[i] = float(line.split(' ')[2])
				self.modes[i,:] = map(float,line.split(' ')[3:])
				self.scales2[i] = float(self.rmsd.get())*np.sqrt(self.atomNumber)/self.scales[i]

			self.arrow_direction = 1

			self._buildMolecule()

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
			w7_1 = Tkinter.Button(group2, text="Make", command=self._drawMode)
			w7_1.grid(row=1,column=2,columnspan=2)
			w7_2 = Tkinter.Button(group2, text="Play", command=self._cleanMode)
			w7_2.grid(row=1,column=4,columnspan=2)
			w7_3 = Tkinter.Button(group2, text="Hide", command=self._hideMode)
			w7_3.grid(row=1,column=6,columnspan=2)
			w7_4 = Tkinter.Button(group2, text="Options", command=self._optionsMode)
			w7_4.grid(row=1,column=8,columnspan=2)

			w8_0 = Tkinter.Label(group2, text="Figures:")
			w8_0.grid(row=2,column=0,columnspan=2)
			w8_1 = Tkinter.Button(group2, text="Clear", command=self._drawMode)
			w8_1.grid(row=2,column=2,columnspan=2)
			w8_2 = Tkinter.Button(group2, text="Close", command=self._cleanMode)
			w8_2.grid(row=2,column=4,columnspan=2)
			w8_3 = Tkinter.Button(group2, text="Hide", command=self._hideMode)
			w8_3.grid(row=2,column=6,columnspan=2)
			w8_4 = Tkinter.Button(group2, text="Options", command=self._optionsMode)
			w8_4.grid(row=2,column=8,columnspan=2)

			w9_0 = Tkinter.Label(group2, text="Molecule:") ## Add the molecule number and check from nmwiz.tcl
			w9_0.grid(row=3,column=0,columnspan=2)
			w9_1 = Tkinter.Button(group2, text="Update", command=self._drawMode)
			w9_1.grid(row=3,column=2,columnspan=2)
			w9_2 = Tkinter.Button(group2, text="Focus", command=self._cleanMode)
			w9_2.grid(row=3,column=4,columnspan=2)
			w9_3 = Tkinter.Button(group2, text="Hide", command=self._hideMode)
			w9_3.grid(row=3,column=6,columnspan=2)
			w9_4 = Tkinter.Button(group2, text="Options", command=self._optionsMode)
			w9_4.grid(row=3,column=8,columnspan=2)

			self._buildArrows()
			self._drawArrows()

			self.molecules = openModels.list()
			self.proteinMol = self.molecules[-2]
			self.arrowMol = self.molecules[-1]
			bonds = self.proteinMol.bonds
			residues = self.proteinMol.residues
			atoms = self.proteinMol.atoms
			for b in bonds:
				b.radius = .2
			for r in residues:
				r.ribbonDisplay = False
			for at in atoms:
				at.display = True

			self._colorArrows()

	def _buildMolecule(self):

		m = Molecule()
		for i in range(self.atomNumber):
			residues = []
			if i == 0:
				r = m.newResidue("ALA"," ",i," ")
				atomCA = m.newAtom("CA",Element("C"))
				atomCA.setCoord(Coord(self.coordinates[i,0],self.coordinates[i,1],self.coordinates[i,2]))
				r.addAtom(atomCA)
				atomCA_old = atomCA
				residues.append(r)
			else:
				r = m.newResidue("ALA"," ",i," ")
				atomCA = m.newAtom("CA",Element("C"))
				atomCA.setCoord(Coord(self.coordinates[i,0],self.coordinates[i,1],self.coordinates[i,2]))
				r.addAtom(atomCA)
				m.newBond(atomCA_old,atomCA)
				atomCA_old = atomCA
				residues.append(r)

		openModels.add([m])

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
		#print self.coordinates[0,:]
		#print self.arrows[0,:]
		self._calcMSF()

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
