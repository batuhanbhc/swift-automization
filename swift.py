#!/usr/bin/python3

# NAME:
#	swift.py
#
# PURPOSE:
#	To automate the processing (post level 1) of the observational data collected from the Swift XRT
#
# EXPLANATION:
#		 
#
# CALLING SEQUENCE:
#
#
# OUTPUTS:
#
#
# EXAMPLES:
#
#
# AUTHORS:
#        K.Zaidi           May, 2018
#		 B.Bahceci

import subprocess, os, re, time
from parameter import *
from pathlib import Path
from astropy.io import fits
import sys
import heasoftpy as hsp

# check the installation #
if not 'HEADAS' in os.environ:
    raise RuntimeError('Heasoft is not initalized')

try:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.abspath(''), '..')))
    import heasoftpy as hsp
except ImportError:
    raise RuntimeError('Please ensure heasoftpy is in your PYTHONPATH')
    
if not hasattr(hsp, 'ftlist'):
    raise RuntimeError('Heasoftpy is not full installed. Please run the install.py script')

cwd = os.getcwd()
script_path = os.path.abspath(__file__)

index = script_path[::-1].find("/")
parentDir = script_path[::-1][index+1:]
parentDir = parentDir[::-1]
##################################################################     PC/WT     ########################################################################################

if ignore_radius >= sourceRadius:
	print("Ignore radius for annulus (inner radius) is bigger than source radius (outer radius)")
	print(f"Enter an ignore radius value that is smaller than {sourceRadius} and rerun the script.")
	quit()

#reads lines of the input.txt file into a list
inputFile = open(parentDir +"/"+ inputTxt,'r')               #open the swift file    
lineNum = 0                                     #initialize the line count
obsdir = []                                 	#initializing the array to hold the location of the observation directories given by swift

for line in inputFile.readlines():              #loop to record the information from the file in 'attributes'
	line = line.strip('\n\r')
	obsdir.append(line)
	lineNum += 1
inputFile.close()

for eachObs in obsdir:
	print("=======================================================================")
	print("=======================================================================")

	ev = '/xrt/event'
	obsdir = eachObs
	evtdir = obsdir + ev

	#creating obsdir_reduced for outdir
	obsdir_temp = obsdir[::-1]
	obsdir_reduced = obsdir_temp[obsdir_temp.find('/'):]
	obsdir_reduced = obsdir_reduced[::-1]

	#creating obsid from obsdir
	obsid_temp = obsdir_temp[:obsdir_temp.find('/')]
	obsid = obsid_temp[::-1]
	print(obsdir)
	
	#output directory for screened files coming out of xrtpipeline
	outdir = cwd + "/"  + obsid + '-xrt'
	print('outdir: ' + outdir)

	#create the outdir directory if it does not already exist
	if Path(outdir).exists() == False:
		os.system("mkdir " + outdir)

	#name of the file where warning messages will be recorded
	warninglog = "warning_" + obsid + ".log"
	warning_file = open(outdir + "/" + warninglog, "w")  
	
	obmode = ""
	PC = False
	WT = False
	#determines if an observation was recorded in a PC or a WT mode
	if (os.path.exists(obsdir)):
		event = os.listdir(evtdir)     		#array holds the foldernames in the observation directory
		
		for file in event:
			if "xwtw2po_cl.evt.gz" in file:
				print('Found observation mode: WT')
				windowed_file = file
				WT = True
			elif "xpcw3po_cl.evt.gz" in file:
				print('Found observation mode: PC')
				photon_file = file
				PC = True

		if PC and WT:
			print("Found both WT and PC mode files.")
			hdu = fits.open(evtdir + "/" + photon_file)
			photon_ontime = hdu[0].header["ONTIME"]
			hdu.close()
			
			hdu = fits.open(evtdir + "/" + windowed_file)
			windowed_ontime = hdu[0].header["ONTIME"]
			hdu.close()

			if windowed_ontime > photon_ontime:
				obmode = "wt"
				print("Current observation mode: WT")
				warning_file.write("Warning: WT mode has been automatically selected over PC mode due to having a larger ontime.\n")
				warning_file.write("WT mode ontime: " + str(windowed_ontime) + " | PC mode ontime: " + str(photon_ontime) + "\n\n")
			else:
				obmode = "pc"
				print("Current observation mode: PC")
				warning_file.write("Warning: PC mode has been automatically selected over WT mode due to having a larger ontime.\n")
				warning_file.write("WT mode ontime: " + str(windowed_ontime) + " | PC mode ontime: " + str(photon_ontime) + "\n\n")

		elif PC:
			print("Observation mode: PC")
			obmode = "pc"
		elif WT:
			print("Observation mode: WT")
			obmode = "wt"
		else:
			print("Neither PC nor WT mode files could be found. Skipping the current observation...")
			continue

	else:
		print('Observation not found')
		continue

	#Find the appropriate event file
	if obmode == "wt":
		xrt_evtfile = "sw" + obsid + "xwtw2po_cl.evt"
	else:
		xrt_evtfile = "sw" + obsid + "xpcw3po_cl.evt"
	
	####################################################################     XRTPIPELINE     ################################################################################

	#recording ObsID
	#obsid = evtdir_split[4] 
	print('ObsID: ' + obsid)

	#Input directory for running xrtpipeline
	indir = obsdir
	print('indir:' + indir)

	#steminputs for xrtpipeline
	steminputs = 'sw' + obsid
	print('steminputs: ' + steminputs)

	#name of the file where the log of xrtpipeline is recorded
	xrtlog = 'xrt_' + obsid + '.log'

	if RA == "":
		srcra = "OBJECT"
	else:
		RA = str(RA)
		srcra = "'" + RA + "'"

	if DEC == "":
		srcdec = "OBJECT"
	else:
		DEC = str(DEC)
		srcdec = "'" + DEC + "'"
	
	xrt = "xrtpipeline indir=" + indir + " steminputs=" + steminputs + " outdir=" + outdir + " srcra=" + srcra + " srcdec=" + srcdec + " createexpomap=yes useexpomap=yes plotdevice=ps correctlc=yes clobber=yes cleanup=no > " + xrtlog
	print('Running pipeline command: ' + xrt)

	#change directory to outdir
	os.chdir(outdir)

	os.system(xrt)
	############################	Finding rmf file	########################################
	if obmode == "wt":
		modeKeyword = "windowed"
		grade = "G0:2"
	elif obmode == "pc":
		modeKeyword = "photon"
		grade = "G0:12"
	else:
		print("Unknown observation mode.\n")
		warning_file.write("Unknown observation mode.\n")
		warning_file.close()
		continue

	print("Running quzcif to find rmf file.\n")
	p = subprocess.Popen( "quzcif SWIFT XRT - - matrix - - datamode.eq."+ modeKeyword +".and.grade.eq." + grade +".and.XRTVSUB.eq.6", stdout=subprocess.PIPE, shell=True)
	(quzcif_output, err) = p.communicate()
	p_status = p.wait()

	quzcif = str(quzcif_output)
	quzcif = quzcif.strip("1b\\n' ").split(" ")
	quzcifFile = quzcif[-1].strip("1b\\n'")
	os.system('cp '+  quzcifFile + ' ' + outdir)
	try:
		print('rmf file location: ' + quzcifFile)
	except Exception as e:
		print(f"Could not find rmf file in CALDB: {e}\n")
		warning_file.write(f"Could not find rmf file in CALDB: {e}\n")
		warning_file.close()
		continue

	xrt_filelist = os.listdir(outdir)

	#finding the rmf file in the xrt files (was copied here from CALDB)
	for rmf in xrt_filelist:
		if (re.search('rmf', rmf, re.M|re.I)):
			rmffile = rmf
	try:
		print('rmf file: ' + rmffile)
	except Exception as e:
		print(f"Could not find rmf file in output directory: {e}\n")
		warning_file.write(f"Could not find rmf file in output directory: {e}\n")
		warning_file.close()
		continue

	########################	Extracting source coordinates	########################
	#finding the region file in the screened xrt files
	for reg in xrt_filelist:
		if ("reg" in reg) and ("po" in reg):
			regfile = reg
			break
	
	#reading the source coordinates from the region file created by the xrtpipeline
	try:
		regfile_cood = open(regfile, 'r')
	except Exception as e:
		print(f"Could not open region file: {e}")
		warning_file.write(f"Could not open region file: {e}")
		warning_file.close()
		continue

	cood = regfile_cood.read()
	regfile_cood.close()
	cood = cood[cood.find('(')+1 : cood.find(')')]
	cood = cood.split(',')
	srcx = cood[0]
	srcx = float(srcx)
	srcy = cood[1]
	srcy = float(srcy)

	#Bin the events in cleaned event file into pixels
	if recalculate_source_center:
		if obmode == "wt":
			try:
				hdu = fits.open("sw" + obsid + "xwtw2po_cl.evt")
			except Exception as e:
				print(f"Exception occured while trying to open sw{obsid}xwtw2po_cl.evt: {e}")
				warning_file.write(f"Exception occured while trying to open sw{obsid}xwtw2po_cl.evt: {e}")
				warning_file.close()
				continue
		
			pixelCounts = {}
			events = hdu[1].data

			for eachEvent in events:
				imageX = eachEvent[1]
				imageY = eachEvent[2]

				coordinate = (imageX, imageY)

				if coordinate in pixelCounts:
					pixelCounts[coordinate] += 1
				else:
					pixelCounts[coordinate] = 1
			
			tempMax = 0
			sourceCoords = (0, 0)
			for coordinate, count in pixelCounts.items():
				if count > tempMax:
					tempMax = count
					sourceCoords = coordinate

			srcx = sourceCoords[0]
			srcy = sourceCoords[1]

			with open(regfile, "w") as tempFile:
				tempFile.write("CIRCLE ("+ str(srcx) +","+ str(srcy) +","+ str(sourceRadius) +")\n")

		elif obmode == "pc":
			with fits.open("sw"+obsid+"xpcw3po_cl.evt") as hdu:
				srcra = hdu[1].header["RA_OBJ"]
				srcdec = hdu[1].header["DEC_OBJ"]

			if Path("position.txt").exists():
				os.system("rm position.txt")
			
			boxHalfWidth = 10	# 10 arcmin
			xrtcentroid = "xrtcentroid infile=sw" + obsid+ "xpcw3po_cl.evt outfile=position.txt outdir=" + outdir + " calcpos=yes interactive=no boxra=" + str(srcra) + " boxdec=" + str(srcdec) + " boxradius=" + str(boxHalfWidth)

			#Running xrtcentroid to refind source position
			os.system(xrtcentroid)

			#Rerun xrtcentroid using previous call's coordinates as input with smaller boxes each iteration
			boxWidths = [5, 2, 1]
			for boxHalfWidth in boxWidths:
				with open("position.txt", "r") as posfile:
					fileLines = posfile.readlines()
					srcra = float(fileLines[2].split(" ")[-1])
					srcdec = float(fileLines[3].split(" ")[-1])
				
				if Path("position.txt").exists():
					os.system("rm position.txt")

				xrtcentroid = "xrtcentroid infile=sw" + obsid+ "xpcw3po_cl.evt outfile=position.txt outdir=" + outdir + " calcpos=yes interactive=no boxra=" + str(srcra) + " boxdec=" + str(srcdec) + " boxradius=" + str(boxHalfWidth)

				#Running xrtcentroid to refind source position
				os.system(xrtcentroid)
		
			with open("position.txt", "r") as posfile:
				fileLines = posfile.readlines()
				srcx = float(fileLines[11].split(" ")[-1])
				srcy = float(fileLines[12].split(" ")[-1])
				source_ra = float(fileLines[2].split(" ")[-1])
				source_dec = float(fileLines[3].split(" ")[-1])

			with open("sw" + obsid + "xpcw3po.reg", "w") as tempFile:
				tempFile.write("CIRCLE ("+ str(srcx) +","+ str(srcy) +","+ str(sourceRadius) +")\n")

	if pileup:
		print("Pileup mode is set to True. Converting source region from circle to annulus.")
		print(f"Inner radius: {ignore_radius}")
		print(f"Outer radius: {sourceRadius}")
		temp_file = open(regfile, "w")
		temp_file.write('annulus (' + str(srcx) + ',' + str(srcy) + ','+str(ignore_radius)+','+str(sourceRadius)+')\n')
		temp_file.close()

	############################################################################################

	# Region manipulation for PC
	if obmode == 'pc':
		srcxr = srcx + 100
		srcxl = srcx - 100
		srcyu = srcy + 100
		srcyd = srcy - 100
		backr = 'sw' + obsid + 'back_right.reg'
		back_right = open(backr, 'w')
		back_right.write('CIRCLE (' + str(srcxr) + ',' + str(srcy) + ','+ str(sourceRadius) +')')
		back_right.close()
		backl = 'sw' + obsid + 'back_left.reg'
		back_left = open(backl, 'w')
		back_left.write('CIRCLE (' + str(srcxl) + ',' + str(srcy) + ','+ str(sourceRadius) +')')
		back_left.close()
		backu = 'sw' + obsid + 'back_up.reg'
		back_up = open(backu, 'w')
		back_up.write('CIRCLE (' + str(srcx) + ',' + str(srcyu) + ','+ str(sourceRadius) +')')
		back_up.close()
		backd = 'sw' + obsid + 'back_down.reg'
		back_down = open(backd, 'w')
		back_down.write('CIRCLE (' + str(srcx) + ',' + str(srcyd) + ','+ str(sourceRadius) + ')')
		back_down.close()
		xsel_PC = 'xsel' + obsid + '_PCback.xco'
		xsel_pcback = open(xsel_PC, 'w')
		xsel_pcback.write('xsel' + obsid + '\nno\n')
		xsel_pcback.write('read e'+ '\n')
		xsel_pcback.write(outdir + '\n')
		xsel_pcback.write(xrt_evtfile + '\n' + 'yes' + '\n')
		xsel_pcback.write('filter region ' + backr + '\nextract spectrum\nyes\n')
		xsel_pcback.write('clear region\nfilter region ' + backl + '\nextract spectrum\nyes\n')
		xsel_pcback.write('clear region\nfilter region ' + backu + '\nextract spectrum\nyes\n')
		xsel_pcback.write('clear region\nfilter region ' + backd + '\nextract spectrum\nyes\n')
		xsel_pcback.write('$cd\n')
		xsel_pcback.write('no\nquit\nyes')
		xsel_pcback.close()
		os.system('xselect @' + xsel_PC)
		back_xpec = []
		counts = []
		back_pcreg = open('xselect.log','r')
		for line in back_pcreg.readlines():
			if (re.search('Spectrum         has', line, re.M|re.I)):
				back_xpec.append(line)
		back_pcreg.close()
		for item in back_xpec:
			item = item[25:item.find('coun')]
			item = str(item).strip()
		item = int(item)
		counts.append(item)
		min_ = min(counts)
		pc_backindex = counts.index(min_)
	# Region manipulation for WT
	elif obmode == 'wt':
		back = 'sw' + obsid + 'back.reg'
		backfile = open(back, 'w')
		backfile.write('annulus (' + str(srcx) + ',' + str(srcy) + ','+str(innerRadius)+','+str(outerRadius)+')\n')
		backfile.close()

	#Writing the .xco file for XSELECT to read the commands from
	print('Event file directory:' + outdir)	#see outdir above

	xsel_filename = 'xsel' + obsid + '.xco'
	xsel_file = open(xsel_filename, 'w')
	xsel_file.write('xsel' + obsid + '\nno\n')
	xsel_file.write('read e'+ '\n')
	xsel_file.write(outdir + '\n')

	#adding further things to the .xco file
	xsel_file.write(xrt_evtfile + '\n' + 'yes' + '\n')
	xsel_file.write('filter region ' + regfile + '\nextract spectrum\nsave spectrum sw' + obsid + '_spectrum.pha\nyes\nclear region\n')
	if obmode == 'pc':
		if pc_backindex == 0:
			xsel_file.write('filter region ' + backr + '\nextract spectrum\nsave spectrum sw' + obsid + 'back_spectrum.pha\nyes\n')
		elif pc_backindex == 1:	
			xsel_file.write('clear region\nfilter region ' + backl + '\nextract spectrum\nsave spectrum sw' + obsid + 'back_spectrum.pha\nyes\n')
		elif pc_backindex == 2:	
			xsel_file.write('clear region\nfilter region ' + backu + '\nextract spectrum\nsave spectrum sw' + obsid + 'back_spectrum.pha\nyes\n')
		elif pc_backindex == 3:
			xsel_file.write('clear region\nfilter region ' + backd + '\nextract spectrum\nsave spectrum sw' + obsid + 'back_spectrum.pha\nyes\n')
	elif obmode == 'wt':
		xsel_file.write('filter region ' + back + '\nextract spectrum\nsave spectrum sw' + obsid + 'back_spectrum.pha\nyes\n')
	#xsel_file.write('$cd\n')
	xsel_file.write('no\nquit\nno')
	xsel_file.close()

	os.system('xselect @' + xsel_filename)
	try:
		if obmode == "wt":
			with fits.open("sw"+obsid+"_spectrum.pha", mode='update') as hdu:
				if pileup:
					hdu[1].header["BACKSCAL"] = sourceRadius*2 - ignore_radius*2
				else:
					hdu[1].header["BACKSCAL"] = sourceRadius*2
				hdu.flush()
			
			with fits.open("sw"+obsid+"xwtw2posr.pha", mode='update') as hdu:
				if pileup:
					hdu[1].header["BACKSCAL"] = sourceRadius*2 - ignore_radius*2
				else:
					hdu[1].header["BACKSCAL"] = sourceRadius*2
				hdu.flush()
			
			with fits.open("sw"+obsid+"back_spectrum.pha", mode='update') as hdu:
				hdu[1].header["BACKSCAL"] = outerRadius - innerRadius - 1
				hdu.flush()
	except Exception as e:
		print(f"Exception occured while trying to update sw{obsid}_spectrum.pha file: {e}\n")
		warning_file.write(f"Exception occured while trying to update sw{obsid}_spectrum.pha file: {e}\n")
		warning_file.close()
		continue

	#finding the arf file in the screened xrt files
	for arf in xrt_filelist:
		if (re.search('arf', arf, re.M|re.I)):
			if (re.search('po', arf, re.M|re.I)):
				arffile = arf
	try:
		print('arf file: ' + arffile)
	except Exception as e:
		print(f"Could not find arf file in output directory: {e}\n")
		warning_file.write(f"Could not find arf file in output directory: {e}\n")
		warning_file.close()
		continue

	# Renew the list variable to to include new files
	xrt_filelist = os.listdir(outdir)

	# Find necessary files for xrtmkarf
	for file in xrt_filelist:
		if "po_ex.img" in file:
			exposurefile = file
		elif "sw" + obsid + "_spectrum.pha" == file:
			spectrumfile = file
		elif "posr.arf" in file:
			oldarffile = file

	# Create arf file using xrtmkarf after spectrum file is created
	# oldarffile variable is just used for naming the output file
	try:
		xrtmkarf = "xrtmkarf outfile=" + oldarffile + " rmffile=" + rmffile + " clobber=yes expofile=" + exposurefile + " phafile=" + spectrumfile + f" srcx={srcx} srcy={srcy} psfflag=yes >> " + xrtlog
		os.system(xrtmkarf)
	except Exception as e:
		print(f"Exception occured while running xrtmkarf command: {e}\n")
		warning_file.write(f"Exception occured while running xrtmkarf command: {e}\n")
		warning_file.close()
		continue

	#grppha
	try:
		grp_out = '!sw' + obsid + '_grp.pha'
		backfile = 'sw' + obsid + 'back_spectrum.pha'
		grppha = "grppha infile= '" + "sw" + obsid + "_spectrum.pha'" + " outfile='" + grp_out + "' chatter=0 comm='group min 10 & bad 0-29 & chkey backfile " + backfile + " & chkey ancrfile " + arffile + " & chkey respfile " + rmffile + " & exit'"
		print(grppha)
		os.system(grppha)
	except Exception as e:
		print(f"Exception occured while creating grppha file: {e}\n")
		warning_file.write(f"Exception occured while creating grppha file: {e}\n")
		warning_file.close()
		continue
	
	warning_file.write("Finished with no errors/warnings.\n")
	warning_file.close()

#########################################################################################################################################################################

if Path(script_path[:script_path.rfind("/")] + "/__pycache__").exists():
	os.system("rm -rf " + script_path[:script_path.rfind("/")] + "/__pycache__")