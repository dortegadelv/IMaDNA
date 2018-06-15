import numpy as np
import scipy, sys, math
from scipy import optimize
import random

Nfeval = 1

def callbackF(Xi):
    global Nfeval
    print '{0:4d}   {1: 6.6f}   {2: 6.6f}'.format(Nfeval, Xi[0], ModelLogLikelihoodOnlyDivTime(Xi))
    Nfeval += 1


def ModelLogLikelihoodOnlyDivTime(Pars):
	import sys, numpy, scipy, math
	from scipy.stats import binom
	from CalculateLLReadDataFunctions_v2 import p_V_G_given_C_E_OriginalPre
	import operator as op
	import array
	from array import array
	import subprocess
	### The next variables should be global and accessible to the function
	# InputFile
	# ContaminationRate
	# AncTime
	# PrefixInput
	# InputProportions
	# OutputSFS
	# RunID
	# OutputFile
	### The variables passed to the function are:
	# N1
	# N2
	# NA
	# TDiv
	# M1
	# M2
	# global CurrentPrefix, InputProportions, OutputSFS, RunID, PrefixInput
	print (CurrentPrefix)
	print (InputProportions)
	print (OutputSFS)
	N1 = 10000
	N2 = 10000
	NA = 10000
	TDiv = Pars[0]
	M1 = 0.00004
	M2 = 0.00004
	print "N1 =",str(N1), " N2 =", str(N2), " NA =", str(NA), " TDiv =", str(TDiv), " M1 =", str(M1), " M2 =", str(M2), " AncT =", str(AncTime)
	F_PrefixInput=PrefixInput
	F_InputProportions=InputProportions
	F_OutputSFS=OutputSFS
	F_RunID = str(RunID)
	F_CurrentPrefix = F_PrefixInput + F_RunID
	# sys.exit("All good here")
	### To run the bash process use cmd = subprocess.Popen(['bash', './FirstLLTestGetSFSWithAncestralMaterial.sh'])
	subprocess.call(["bash", "./GetSFSWithAncestralMaterialNoMigrationDifMigPar.sh", str(N1), str(N2), str(NA), str(TDiv), str(M1), str(M2), str(AncTime), F_CurrentPrefix, F_InputProportions, F_OutputSFS])
	# sys.exit("Until here")
	SFSFileName = OutputSFS
	SiteLLFileName = InputFile
	# print ("Contamination rate = " + str(ContaminationRate))
	SFSFile = open (OutputSFS, 'r')
	ThisSFS = numpy.zeros((3,3))
	RowNumber = 0
	for line in SFSFile:
		# print(line)
		SFS = line.split('\t')
		for j in range(0, 3):
			ThisSFS[RowNumber][j] = SFS[j]
		RowNumber = RowNumber + 1
		if (RowNumber == 3):
			break
	File = open (SiteLLFileName, 'r')
	FullLikelihood = 0.0
	RowNumber = 0
	BinomVector = numpy.zeros(20)
	CombinVector = numpy.zeros(5)
	for line in File:
		#print(line)
		LineComponents = line.split('\t')
		PresentDayDerivedAlleleCount = int(LineComponents[0])
		AlleleFrequency = float(LineComponents[1])
		AncestralAlleleCounts = int(LineComponents[2])
		DerivedAlleleCounts = int(LineComponents[3])
		EPlus = float(LineComponents[4])
		EMinus = float(LineComponents[5])
		EPlusPlus = float(LineComponents[6])
		EMinusMinus = float(LineComponents[7])
		EPlusMinus = float(LineComponents[8])
		EFullMatrixXMatrix = float(LineComponents[9])
		ECero = 1.0 - EPlus - EMinus - EFullMatrixXMatrix
		R_i = AncestralAlleleCounts + DerivedAlleleCounts
		u_i = DerivedAlleleCounts
		BinomVector[0] = float(0.0)
		BinomVector[1] = float(0.0)
		BinomVector[2] = float(0.0)
		BinomVector[3] = float(0.0)
		BinomVector[4] = float(0.0)
		BinomVector[5] = float(0.0)
		BinomVector[6] = float(0.0)
		BinomVector[7] = float(0.0)
		BinomVector[8] = float(0.0)
		BinomVector[9] = float(0.0)
		BinomVector[10] = float(0.0)
		BinomVector[11] = float(0.0)
		BinomVector[12] = float(0.0)
		BinomVector[13] = float(0.0)
		BinomVector[14] = float(0.0)
		BinomVector[15] = float(0.0)
		BinomVector[16] = float(0.0)
		BinomVector[17] = float(0.0)
		BinomVector[18] = float(0.0)
		BinomVector[19] = float(0.0)
		CombinVector[0] = float(LineComponents[10])
		CombinVector[1] = float(LineComponents[11])
		CombinVector[2] = float(LineComponents[12])
		CombinVector[3] = float(LineComponents[13])
		CombinVector[4] = float(LineComponents[14])
		BinomVectorA = array('d',BinomVector)
		CombinVectorA = array('d',CombinVector)
		### Full Likelihood
		RowNumber = RowNumber + 1
		# print(type(PresentDayDerivedAlleleCount))
		# print(type(AlleleFrequency))
		# print(type(AncestralAlleleCounts))
		# print(type(DerivedAlleleCounts))
		# print(type(EPlus))
		# print(type(EMinus))
		# print(type(EPlusPlus))
		# print(type(EMinusMinus))
		# print(type(EPlusMinus))
		# print(type(ECero))
		# print(type(ContaminationRate))
		# FullLikelihood = FullLikelihood + math.log(p_ui_Gi_Original( DerivedAlleleCounts, AncestralAlleleCounts, PresentDayDerivedAlleleCount, ThisSFS))
		# FullLikelihood = FullLikelihood + math.log(p_ui_Gi_C_Original( PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts, DerivedAlleleCounts, ThisSFS, ContaminationRate))
		FullLikelihood = FullLikelihood + math.log(p_V_G_given_C_E_OriginalPre(PresentDayDerivedAlleleCount,AlleleFrequency,AncestralAlleleCounts,DerivedAlleleCounts,EPlus,EMinus,EPlusPlus,EMinusMinus,EPlusMinus,ECero,ThisSFS,ContaminationRate,BinomVectorA,CombinVectorA))
		if (RowNumber % 10000000 == 0):
			print (str(RowNumber))
	return(- FullLikelihood)


def ModelLogLikelihood(Pars):
	import sys, numpy, scipy, math
	from scipy.stats import binom
	from CalculateLLReadDataFunctions_v2 import p_V_G_given_C_E_OriginalPre
	import operator as op
	import array
	from array import array
	import subprocess
	### The next variables should be global and accessible to the function
	# InputFile
	# ContaminationRate
	# AncTime
	# PrefixInput
	# InputProportions
	# OutputSFS
	# RunID
	# OutputFile
	### The variables passed to the function are:
	# N1
	# N2
	# NA
	# TDiv
	# M1
	# M2
	# global CurrentPrefix, InputProportions, OutputSFS, RunID, PrefixInput
	print (CurrentPrefix)
	print (InputProportions)
	print (OutputSFS)
	N1 = round(Pars[0])
	N2 = round(Pars[1])
	NA = round(Pars[2])
	TDiv = round(Pars[3])
	M1 = round(Pars[4],8)
	M2 = round(Pars[5],8)
	print "N1 =",str(N1), " N2 =", str(N2), " NA =", str(NA), " TDiv =", str(TDiv), " M1 =", str(M1), " M2 =", str(M2), " AncT =", str(AncTime)
	F_PrefixInput=PrefixInput
	F_InputProportions=InputProportions
	F_OutputSFS=OutputSFS
	F_RunID = str(RunID)
	F_CurrentPrefix = F_PrefixInput + F_RunID
	# sys.exit("All good here")
	### To run the bash process use cmd = subprocess.Popen(['bash', './FirstLLTestGetSFSWithAncestralMaterial.sh'])
	subprocess.call(["bash", "./GetSFSWithAncestralMaterialNoMigration.sh", str(N1), str(N2), str(NA), str(TDiv), str(M1), str(M2), str(AncTime), F_CurrentPrefix, F_InputProportions, F_OutputSFS])
	# sys.exit("Until here")
	SFSFileName = OutputSFS
	SiteLLFileName = InputFile
	print ("Contamination rate = " + str(ContaminationRate))
	# SFSFile = open (OutputSFS, 'r')
	# ThisSFS = numpy.zeros((3,3))
	# RowNumber = 0
	# for line in SFSFile:
		# print(line)
		# SFS = line.split('\t')
		# for j in range(0, 3):
			# ThisSFS[RowNumber][j] = SFS[j]
		# RowNumber = RowNumber + 1
		# if (RowNumber == 3):
			# break
	File = open (SiteLLFileName, 'r')
	FullLikelihood = 0.0
	RowNumber = 0
	BinomVector = numpy.zeros(20)
	CombinVector = numpy.zeros(5)
	for line in File:
		#print(line)
		LineComponents = line.split('\t')
		PresentDayDerivedAlleleCount = int(LineComponents[0])
		AlleleFrequency = float(LineComponents[1])
		AncestralAlleleCounts = int(LineComponents[2])
		DerivedAlleleCounts = int(LineComponents[3])
		EPlus = float(LineComponents[4])
		EMinus = float(LineComponents[5])
		EPlusPlus = float(LineComponents[6])
		EMinusMinus = float(LineComponents[7])
		EPlusMinus = float(LineComponents[8])
		EFullMatrixXMatrix = float(LineComponents[9])
		ECero = 1.0 - EPlus - EMinus - EFullMatrixXMatrix
		R_i = AncestralAlleleCounts + DerivedAlleleCounts
		u_i = DerivedAlleleCounts
		BinomVector[0] = float(0.0)
		BinomVector[1] = float(0.0)
		BinomVector[2] = float(0.0)
		BinomVector[3] = float(0.0)
		BinomVector[4] = float(0.0)
		BinomVector[5] = float(0.0)
		BinomVector[6] = float(0.0)
		BinomVector[7] = float(0.0)
		BinomVector[8] = float(0.0)
		BinomVector[9] = float(0.0)
		BinomVector[10] = float(0.0)
		BinomVector[11] = float(0.0)
		BinomVector[12] = float(0.0)
		BinomVector[13] = float(0.0)
		BinomVector[14] = float(0.0)
		BinomVector[15] = float(0.0)
		BinomVector[16] = float(0.0)
		BinomVector[17] = float(0.0)
		BinomVector[18] = float(0.0)
		BinomVector[19] = float(0.0)
		CombinVector[0] = float(LineComponents[10])
		CombinVector[1] = float(LineComponents[11])
		CombinVector[2] = float(LineComponents[12])
		CombinVector[3] = float(LineComponents[13])
		CombinVector[4] = float(LineComponents[14])
		BinomVectorA = array('d',BinomVector)
		CombinVectorA = array('d',CombinVector)
		### Full Likelihood
		RowNumber = RowNumber + 1
		# print(type(PresentDayDerivedAlleleCount))
		# print(type(AlleleFrequency))
		# print(type(AncestralAlleleCounts))
		# print(type(DerivedAlleleCounts))
		# print(type(EPlus))
		# print(type(EMinus))
		# print(type(EPlusPlus))
		# print(type(EMinusMinus))
		# print(type(EPlusMinus))
		# print(type(ECero))
		# print(type(ContaminationRate))
		# FullLikelihood = FullLikelihood + math.log(p_ui_Gi_Original( DerivedAlleleCounts, AncestralAlleleCounts, PresentDayDerivedAlleleCount, ThisSFS))
		# FullLikelihood = FullLikelihood + math.log(p_ui_Gi_C_Original( PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts, DerivedAlleleCounts, ThisSFS, ContaminationRate))
		FullLikelihood = FullLikelihood + math.log(p_V_G_given_C_E_OriginalPre(PresentDayDerivedAlleleCount,AlleleFrequency,AncestralAlleleCounts,DerivedAlleleCounts,EPlus,EMinus,EPlusPlus,EMinusMinus,EPlusMinus,ECero,ThisSFS,ContaminationRate,BinomVectorA,CombinVectorA))
		if (RowNumber % 10000000 == 0):
			print (str(RowNumber))
	return(-FullLikelihood)

def ModelLogLikelihoodDifParameterization(Pars):
	import sys, numpy, scipy, math
	from scipy.stats import binom
	from CalculateLLReadDataFunctions_v2 import p_V_G_given_C_E_OriginalPreDataRead
	import operator as op
	import array
	from array import array
	import subprocess
	### The next variables should be global and accessible to the function
	# InputFile
	# ContaminationRate
	# AncTime
	# PrefixInput
	# InputProportions
	# OutputSFS
	# RunID
	# OutputFile
	### The variables passed to the function are:
	# N1
	# N2
	# NA
	# TDiv
	# M1
	# M2
	# global CurrentPrefix, InputProportions, OutputSFS, RunID, PrefixInput
	print (CurrentPrefix)
	print (InputProportions)
	print (OutputSFS)
	N1 = 10** Pars[0]
	N2 = 10 ** Pars[1]
	NA = 10 ** Pars[2]
	TDiv = 10 ** Pars[3]
	M1 = 10 ** Pars[4]
	M2 = 10 ** Pars[5]
	print "N1 =",str(N1), " N2 =", str(N2), " NA =", str(NA), " TDiv =", str(TDiv), " M1 =", str(M1), " M2 =", str(M2), " AncT =", str(AncTime) , " Contamination Rate = ", str(ContaminationRate)
	F_PrefixInput=PrefixInput
	F_InputProportions=InputProportions
	F_OutputSFS=OutputSFS
	F_RunID = str(RunID)
	F_CurrentPrefix = F_PrefixInput + F_RunID
	F_CoalDivNumber = CoalDivNumber
	# print "Coal Div Number = " + str(CoalDivNumber)
	# sys.exit("All good here")
	### To run the bash process use cmd = subprocess.Popen(['bash', './FirstLLTestGetSFSWithAncestralMaterial.sh'])
	subprocess.call(["bash", "/opt/IMaDNA/GetSFSWithAncestralMaterialNoMigration.sh", str(N1), str(N2), str(NA), str(TDiv), str(M1), str(M2), str(AncTime), F_CurrentPrefix, F_InputProportions, F_OutputSFS, str(F_CoalDivNumber)])
	print("bash /opt/IMaDNA/GetSFSWithAncestralMaterialNoMigration.sh" + " " + str(N1) + " " + str(N2) + " " + str(NA) + " " + str(TDiv) + " " + str(M1) + " " + str(M2) + " " + str(AncTime) + " " + str(F_CurrentPrefix) + " " + str(F_InputProportions) + " " + str(F_OutputSFS) + " " + str(F_CoalDivNumber))
	# sys.exit("Until here")
	SFSFileName = OutputSFS
	# SiteLLFileName = InputFile
	# print ("Contamination rate = " + str(ContaminationRate))
	SFSFile = open (OutputSFS, 'r')
	ThisSFS = numpy.zeros((3,3))
	RowNumber = 0
	for line in SFSFile:
		print(line)
		SFS = line.split('\t')
		for j in range(0, 3):
			ThisSFS[RowNumber][j] = SFS[j]
		RowNumber = RowNumber + 1
		if (RowNumber == 3):
			break
	# File = open (SiteLLFileName, 'r')
	FullLikelihood = 0.0
	RowNumber = 0
	BinomVector = numpy.zeros(20)
	CombinVector = numpy.zeros(5)
	StartingRow = 0
	print "Number of rows = " + str(len(SitesInformation)) + "\n"
	for line in range(0,len(SitesInformation)):
		#print(line)
		PresentDayDerivedAlleleCount = int(SitesInformation[line][1])
		AlleleFrequency = float(SitesInformation[line][2])
		AncestralAlleleCounts = int(SitesInformation[line][3])
		DerivedAlleleCounts = int(SitesInformation[line][4])
		# EPlus = float(LineComponents[5])
		# EMinus = float(LineComponents[6])
		# EPlusPlus = float(LineComponents[7])
		# EMinusMinus = float(LineComponents[8])
		# EPlusMinus = float(LineComponents[9])
		# EFullMatrixXMatrix = float(LineComponents[10])
		# ECero = 1.0 - EPlus - EMinus - EFullMatrixXMatrix
		R_i = AncestralAlleleCounts + DerivedAlleleCounts
		u_i = DerivedAlleleCounts
		# BinomVector = []
		SNPTypeNumber = int(SitesInformation[line][0])
		# print "Line = " + str(RowNumber) + "\t" + str(len(LineComponents))
		for Element in range(10, 25 ):
			# print "Element = " + str(Element) + " Value = " + str(LineComponents[Element])
			BinomVector[Element - 10] = float(SitesInformation[line][Element])
		# print "Line = " + str(RowNumber) + "\t" + str(len(LineComponents))
		CombinVector[0] = float(SitesInformation[line][5])
		CombinVector[1] = float(SitesInformation[line][6])
		CombinVector[2] = float(SitesInformation[line][7])
		CombinVector[3] = float(SitesInformation[line][8])
		CombinVector[4] = float(SitesInformation[line][9])
		BinomVectorA = array('d',BinomVector)
		CombinVectorA = array('d',CombinVector)
		### Full Likelihood
		RowNumber = RowNumber + 1
		# print(type(PresentDayDerivedAlleleCount))
		# print(type(AlleleFrequency))
		# print(type(AncestralAlleleCounts))
		# print(type(DerivedAlleleCounts))
		# print(type(EPlus))
		# print(type(EMinus))
		# print(type(EPlusPlus))
		# print(type(EMinusMinus))
		# print(type(EPlusMinus))
		# print(type(ECero))
		# print(type(ContaminationRate))
		# FullLikelihood = FullLikelihood + math.log(p_ui_Gi_Original( DerivedAlleleCounts, AncestralAlleleCounts, PresentDayDerivedAlleleCount, ThisSFS))
		# FullLikelihood = FullLikelihood + math.log(p_ui_Gi_C_Original( PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts, DerivedAlleleCounts, ThisSFS, ContaminationRate))
		CurrentProbabilities = p_V_G_given_C_E_OriginalPreDataRead(PresentDayDerivedAlleleCount,AlleleFrequency,AncestralAlleleCounts,DerivedAlleleCounts,ThisSFS,ContaminationRate,BinomVectorA,CombinVectorA)
		# print str(CurrentProbabilities[0]) + "\t" + str(CurrentProbabilities[1]) + "\t" + str(CurrentProbabilities[2]) + "\t" + str(CurrentProbabilities[3]) + "\t" + str(CurrentProbabilities[4]) + "\n"
		for ErrorElements in range(StartingRow,StartingRow + SNPTypeNumber):
			ThisProb = 0
			ECero = 1.0 - ErrorInformation[ErrorElements][0] - ErrorInformation[ErrorElements][1] - ErrorInformation[ErrorElements][5]
			ThisProb = (ErrorInformation[ErrorElements][4] + ECero) * CurrentProbabilities[0]
			ThisProb = ThisProb + ErrorInformation[ErrorElements][3] * CurrentProbabilities[1]
			ThisProb = ThisProb + ErrorInformation[ErrorElements][1] * CurrentProbabilities[2]
			ThisProb = ThisProb + ErrorInformation[ErrorElements][0] * CurrentProbabilities[3]
			ThisProb = ThisProb + ErrorInformation[ErrorElements][2] * CurrentProbabilities[4]
			if (ThisProb != 0.0):
				FullLikelihood = FullLikelihood + math.log(ThisProb)
		# FullLikelihood = FullLikelihood + SNPTypeNumber * math.log(p_V_G_given_C_E_OriginalPre(PresentDayDerivedAlleleCount,AlleleFrequency,AncestralAlleleCounts,DerivedAlleleCounts,EPlus,EMinus,EPlusPlus,EMinusMinus,EPlusMinus,ECero,ThisSFS,ContaminationRate,BinomVectorA,CombinVectorA))
		StartingRow = StartingRow + SNPTypeNumber
		if (RowNumber % 100000000 == 0):
			print (str(RowNumber))
	return(-FullLikelihood)

CurrentPrefix='Prefix'
PrefixInput='Prefix'
InputProportions='OutputProportions.txt'
OutputSFS='OutputSFS.txt'
ParameterFile = open (sys.argv[1],'r')
RunID=sys.argv[2]
N1_startflag = 0
N2_startflag = 0
NA_startflag = 0
DivTime_startflag = 0
FourNm1_startflag = 0
FourNm2_startflag = 0
LowN1Boundary = 100
UpN1Boundary = 100000
LowN2Boundary = 100
UpN2Boundary = 100000
LowNABoundary = 100
UpNABoundary = 100000
LowDivTimeBoundary = 10030
UpDivTimeBoundary = 100000
LowFourNm1Boundary = 0.0002
UpFourNm1Boundary = 10
LowFourNm2Boundary = 0.0002
UpFourNm2Boundary = 10
Iterations = 1
AncTime=10000
CoalDivNumber = 1000
AllIterationsFile = 'PrintAllIterations.txt'
AllIterationsFlag = 0
HeaderFile = "N1_initial\tN2_initial\tNA_initial\tDivTime_initial\t4Nm1_initial\t4Nm2_initial\tN1_final\tN2_final\tNA_final\tDivTime_final\t4Nm1_final\t4Nm2_final\tLog-likelihood\n"

for line in ParameterFile:
	Test = line.rstrip()
	if len(Test) > 0:
		Test2 = Test.split()
		Test2Lower = Test2[0].lower()
		if (Test2Lower == 'errorfile'):
			ErrorFile=Test2[1]
		if (Test2Lower == 'inputfile'):
			InputFile=Test2[1]
		if (Test2Lower == 'contaminationrate'):
			ContaminationRate=float(Test2[1])
		if (Test2Lower == 'anctime'):
			AncTime=int(Test2[1])
		if (Test2Lower == 'outputresults'):
			OutputResults=Test2[1]
		if (Test2Lower == 'prefix'):
			CurrentPrefix = Test2[1]
			PrefixInput = Test2[1]
		if (Test2Lower == 'inputproportions'):
			InputProportions = Test2[1]
		if (Test2Lower == 'sfs'):
			OutputSFS= Test2[1]
		if (Test2Lower == 'runid'):
			RunID = Test2[1]
		if (Test2Lower == 'startn1'):
			N1_start = Test2[1]
			N1_startflag = 1
		if (Test2Lower == 'startn2'):
			N2_start = Test2[1]
			N2_startflag = 1
		if (Test2Lower == 'startna'):
			NA_start = Test2[1]
			NA_startflag = 1
		if (Test2Lower == 'startdivtime'):
			DivTime_start = Test2[1]
			DivTime_startflag = 1
		if (Test2Lower == 'start4nm1'):
			FourNm1_start = Test2[1]
			FourNm1_startflag = 1
		if (Test2Lower == 'start4nm2'):
			FourNm2_start = Test2[1]
			FourNm2_startflag = 1
		if (Test2Lower == 'n1boundaries'):
			LowN1Boundary = float(Test2[1])
			UpN1Boundary = float(Test2[2])
			if (LowN1Boundary > UpN1Boundary):
				print "Error Low N1 boundary " + str(LowN1Boundary) + " is bigger than upper boundary" + str(UpN1Boundary)
				sys.exit()
		if (Test2Lower == 'n2boundaries'):
			LowN2Boundary = float(Test2[1])
			UpN2Boundary = float(Test2[2])
			if (LowN2Boundary > UpN2Boundary):
				print "Error Low N2 boundary " + str(LowN2Boundary) + " is bigger than upper boundary" + str(UpN2Boundary)
				sys.exit()
		if (Test2Lower == 'naboundaries'):
			LowNABoundary = float(Test2[1])
			UpNABoundary = float(Test2[2])
			if (LowNABoundary > UpNABoundary):
				print "Error Low NA boundary " + str(LowNABoundary) + " is bigger than upper boundary" + str(UpNABoundary)
				sys.exit()
		if (Test2Lower == 'divtimeboundaries'):
			LowDivTimeBoundary = float(Test2[1])
			UpDivTimeBoundary = float(Test2[2])
			if (LowDivTimeBoundary > UpDivTimeBoundary):
				print "Error Low DivTime boundary " + str(LowDivTimeBoundary) + " is bigger than upper boundary" + str(UpDivTimeBoundary)
				sys.exit()
		if (Test2Lower == '4nm1boundaries'):
			LowFourNm1Boundary = float(Test2[1])
			UpFourNm1Boundary = float(Test2[2])
			if (LowFourNm1Boundary > UpFourNm1Boundary):
				print "Error Low FourNm1 boundary " + str(LowFourNm1Boundary) + " is bigger than upper boundary" + str(UpFourNm1Boundary)
				sys.exit()
		if (Test2Lower == '4nm2boundaries'):
			LowFourNm2Boundary = float(Test2[1])
			UpFourNm2Boundary = float(Test2[2])
			if (LowFourNm2Boundary > UpFourNm2Boundary):
				print "Error Low FourNm2 boundary " + str(LowFourNm2Boundary) + " is bigger than upper boundary" + str(UpFourNm2Boundary)
				sys.exit()
		if (Test2Lower == 'iterations'):
			Iterations = int(Test2[1])
		if (Test2Lower == 'alliterationsfile'):
			AllIterationsFlag = 1
			AllIterationsFile = Test2[1]
			IterationsPrintFile = open (AllIterationsFile, 'w')
			IterationsPrintFile.write(HeaderFile)
			IterationsPrintFile.flush()
OutputFileResults = open (OutputResults, 'w')

##### Initialization

SiteLLFileName = InputFile
File = open (SiteLLFileName, 'r')
SitesInformation = []
for line in File:
	CurrentLine = []
	LineComponents = line.split()
	CurrentLine.append(int(LineComponents[0]))
	CurrentLine.append(int(LineComponents[1]))
	CurrentLine.append(float(LineComponents[2]))
	CurrentLine.append(int(LineComponents[3]))
	CurrentLine.append(int(LineComponents[4]))
	CurrentLine.append(float(LineComponents[5]))
	CurrentLine.append(float(LineComponents[6]))
	CurrentLine.append(float(LineComponents[7]))
	CurrentLine.append(float(LineComponents[8]))
	CurrentLine.append(float(LineComponents[9]))
	for Element in range(10, len(LineComponents) ):
		CurrentLine.append(float(LineComponents[Element]))
	SitesInformation.append(CurrentLine)
File.close()
ErrorInformation = []
File = open (ErrorFile, 'r')
for line in File:
	CurrentLine = []
	LineComponents = line.split()
	CurrentLine.append(float(LineComponents[0]))
	CurrentLine.append(float(LineComponents[1]))
	CurrentLine.append(float(LineComponents[2]))
	CurrentLine.append(float(LineComponents[3]))
	CurrentLine.append(float(LineComponents[4]))
	CurrentLine.append(float(LineComponents[5]))
	ErrorInformation.append(CurrentLine)

for i in range(0, Iterations):
	if (N1_startflag != 1):
		N1_start = random.uniform(math.log(LowN1Boundary,10),math.log(UpN1Boundary,10))
	if (N2_startflag != 1):
		N2_start = random.uniform(math.log(LowN2Boundary,10),math.log(UpN2Boundary,10))
	if (NA_startflag != 1):
		NA_start = random.uniform(math.log(LowNABoundary,10),math.log(UpNABoundary,10))
	if (DivTime_startflag != 1):
		DivTime_start = random.uniform(math.log(LowDivTimeBoundary,10),math.log(UpDivTimeBoundary,10))
	if (FourNm1_startflag != 1):
		FourNm1_start = random.uniform(math.log(LowFourNm1Boundary,10),math.log(UpFourNm1Boundary,10))
	if (FourNm2_startflag != 1):
		FourNm2_start = random.uniform(math.log(LowFourNm2Boundary,10),math.log(UpFourNm2Boundary,10))
	
	InverseLowN1Boundary = 1. / LowN1Boundary
	InverseLowN2Boundary = 1. / LowN2Boundary
	InverseLowNABoundary = 1. / LowNABoundary
	InverseLowDivTimeBoundary = 1. / LowDivTimeBoundary
	
	InverseUpN1Boundary = 1. / UpN1Boundary
	InverseUpN2Boundary = 1. / UpN2Boundary
	InverseUpNABoundary = 1. / UpNABoundary
	InverseUpDivTimeBoundary = 1. / UpDivTimeBoundary
	
	
	x0 = np.array([ N1_start, N2_start, NA_start , DivTime_start , FourNm1_start, FourNm2_start])
	print(x0)
	Data = scipy.optimize.minimize(ModelLogLikelihoodDifParameterization,x0, method = 'L-BFGS-B',bounds =((math.log(LowN1Boundary,10), math.log(UpN1Boundary,10)), (math.log(LowN2Boundary,10), math.log(UpN2Boundary,10)),(math.log(LowNABoundary,10), math.log(UpNABoundary,10)),(math.log(LowDivTimeBoundary,10), math.log(UpDivTimeBoundary,10)),(math.log(LowFourNm1Boundary,10), math.log(UpFourNm1Boundary,10)),(math.log(LowFourNm2Boundary,10), math.log(UpFourNm2Boundary,10)),), options = {'eps' : 1e-4 })
	
	InitialValues = x0.tolist()
	FinalValues = Data.x.tolist()
	LL = [Data.fun]
	Array = InitialValues + FinalValues + LL
	if (AllIterationsFlag == 1):
		for j in range(0,12):
			Value = 10**Array[j]
			IterationsPrintFile.write(str(Value))
			IterationsPrintFile.write("\t")
		IterationsPrintFile.write(str(Array[12]))
		IterationsPrintFile.write("\n")
		IterationsPrintFile.flush()
	if (i == 0):
		BestLL = Data.fun
		BestArray = Array
	else:
		if (Data.fun < BestLL):
			BestLL = Data.fun
			BestArray = Array
	# OutputFileResults.write(str(Array))
if (AllIterationsFlag == 1):
	IterationsPrintFile.close()
OutputFileResults.write(HeaderFile) 
for j in range(0,12):
	Value = 10**BestArray[j]
	OutputFileResults.write(str(Value))
	OutputFileResults.write("\t")
OutputFileResults.write(str(BestArray[12]))
OutputFileResults.write("\n")
OutputFileResults.close()

