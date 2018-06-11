# cimport numpy
cimport cython
import sys, scipy, math
from scipy.stats import binom
import operator as op
from cpython cimport array
import array

def PrintTwoNumbers (int[:] Numbers):
	# cdef int *my_ints
	# my_ints = <int *>malloc(len(a)*cython.sizeof(int))
	# for i in xrange(len(a)):
	print(Numbers[0])
	print(Numbers[1])

def primes(int kmax):
	cdef int n, k, i
	cdef int p[1000]
	result = []
	if kmax > 1000:
		kmax = 1000
	k = 0
	n = 2
	while k < kmax:
		i = 0
		while i < k and n % p[i] != 0:
			i = i + 1
		if i == k:
			p[k] = n
			k = k + 1
			result.append(n)
		n = n + 1
	return result

def ncr_OriginalPre(n, r):
	r = min(r, n-r)
	if r == 0: return 1
	numer = reduce(op.mul, xrange(n, n-r, -1))
	denom = reduce(op.mul, xrange(1, r+1))
	return numer//denom

def p_V_G_given_C_E_OriginalPre(int PresentDayDerivedAlleleCount, double AlleleFrequency, int AncestralAlleleCounts, int DerivedAlleleCounts, double EPlus, double EMinus, double EPlusPlus, double EMinusMinus, double EPlusMinus, double ECero, double[:,:] ThisSFS, double ContaminationRate, double[:] BinomVector, double[:] CombinVector):
	cdef int R_i = AncestralAlleleCounts + DerivedAlleleCounts
	if R_i == 0:
		return 1
	cdef double ProbabilityPartTwo, ProbabilityPartThree, ProbabilityPartFour, ProbabilityPartFive
	cdef int CurIndex = 0
	cdef double ProbabilityPartOne = (ECero + EPlusMinus) * p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts, DerivedAlleleCounts, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[2]
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1) * (AncestralAlleleCounts + 1)
	CurIndex = CurIndex + 3
	if (CombinVector[0] != 0) and (EMinusMinus != 0.0):
		ProbabilityPartTwo = ( EMinusMinus) * p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts + 2, DerivedAlleleCounts - 2, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[0]
	else:
		ProbabilityPartTwo = 0.0
	CurIndex = CurIndex + 3
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1 - 2) * (AncestralAlleleCounts + 1 + 2)
	if (CombinVector[1] != 0) and (EMinus != 0.0):
		ProbabilityPartThree = (EMinus) * p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts + 1, DerivedAlleleCounts - 1, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[1]
	else:
		ProbabilityPartThree = 0.0
	CurIndex = CurIndex + 3
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1 - 1) * (AncestralAlleleCounts + 1 + 1)
	if (CombinVector[3] != 0) and (EPlus != 0.0):
		ProbabilityPartFour = (EPlus) * p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts - 1, DerivedAlleleCounts + 1, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[3]
	else:
		ProbabilityPartFour = 0.0
	CurIndex = CurIndex + 3
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1 + 1) * (AncestralAlleleCounts + 1 - 1)
	if (CombinVector[4] != 0) and (EPlusPlus != 0.0):
		ProbabilityPartFive = (EPlusPlus) * p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts - 2, DerivedAlleleCounts + 2, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[4]
	else:
		ProbabilityPartFive = 0.0
	cdef double Probability = ProbabilityPartOne + ProbabilityPartTwo + ProbabilityPartThree + ProbabilityPartFour + ProbabilityPartFive
	# print "ProbPartOne = " + str(ProbabilityPartOne) + " ProbPartTwo = " + str(ProbabilityPartTwo) + " ProbPartThree = " + str(ProbabilityPartThree) + " ProbPartFour = " + str(ProbabilityPartFour) + " ProbPartFive = " + str(ProbabilityPartFive) + " ProbTotal = " + str(Probability)
	# print "CombinationsPartOne = " + str(CombinVector[0]) + " CombinationsPartTwo = " + str(CombinVector[1]) + " CombinationsPartThree = " + str(CombinVector[2]) + " CombinationsPartFour = " + str(CombinVector[3]) + " CombinationsPartFive = " + str(CombinVector[4])
	return Probability

def p_V_G_given_C_E_OriginalPreDataRead(int PresentDayDerivedAlleleCount, double AlleleFrequency, int AncestralAlleleCounts, int DerivedAlleleCounts, double[:,:] ThisSFS, double ContaminationRate, double[:] BinomVector, double[:] CombinVector):
	cdef int R_i = AncestralAlleleCounts + DerivedAlleleCounts
	cpdef double Probability[5]
	if R_i == 0:
		Probability[0] = 0.0
		Probability[1] = 0.0
		Probability[2] = 0.0
		Probability[3] = 0.0
		Probability[4] = 0.0
		return Probability
	# print " Anc = " + str(AncestralAlleleCounts) + " Der = " + str(DerivedAlleleCounts) + "Ri = " + str(R_i) + "\n"
	cdef double ProbabilityPartTwo, ProbabilityPartThree, ProbabilityPartFour, ProbabilityPartFive
	cdef int CurIndex = 0
	# cpdef double Probability[5]
	cdef double ProbabilityPartOne = p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts, DerivedAlleleCounts, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[2]
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1) * (AncestralAlleleCounts + 1)
	CurIndex = CurIndex + 3
	if (CombinVector[0] != 0):
		ProbabilityPartTwo = p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts + 2, DerivedAlleleCounts - 2, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[0]
	else:
		ProbabilityPartTwo = 0.0
	CurIndex = CurIndex + 3
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1 - 2) * (AncestralAlleleCounts + 1 + 2)
	if (CombinVector[1] != 0):
		ProbabilityPartThree = p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts + 1, DerivedAlleleCounts - 1, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[1]
	else:
		ProbabilityPartThree = 0.0
	CurIndex = CurIndex + 3
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1 - 1) * (AncestralAlleleCounts + 1 + 1)
	if (CombinVector[3] != 0):
		ProbabilityPartFour = p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts - 1, DerivedAlleleCounts + 1, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[3]
	else:
		ProbabilityPartFour = 0.0
	CurIndex = CurIndex + 3
	# CurIndex = CurIndex + (DerivedAlleleCounts + 1 + 1) * (AncestralAlleleCounts + 1 - 1)
	if (CombinVector[4] != 0):
		ProbabilityPartFive = p_ui_Gi_C_OriginalPre (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts - 2, DerivedAlleleCounts + 2, ThisSFS, ContaminationRate, BinomVector, CurIndex, CombinVector)/ CombinVector[4]
	else:
		ProbabilityPartFive = 0.0
	Probability[0] = ProbabilityPartOne
	Probability[1] = ProbabilityPartTwo
	Probability[2] = ProbabilityPartThree
	Probability[3] = ProbabilityPartFour
	Probability[4] = ProbabilityPartFive
	# print str(Probability[0]) + "\t" + str(Probability[1]) + "\t" + str(Probability[2]) + "\t" + str(Probability[3]) + "\t" + str(Probability[4]) + "\n"
	# print str(ProbabilityPartOne)+ "\t" + str(ProbabilityPartTwo) + "\t" + str(ProbabilityPartThree) + "\t" + str(ProbabilityPartFour) + "\t" + str(ProbabilityPartFive) + "\n"
	# cdef double Probability = ProbabilityPartOne + ProbabilityPartTwo + ProbabilityPartThree + ProbabilityPartFour + ProbabilityPartFive
	# print "ProbPartOne = " + str(ProbabilityPartOne) + " ProbPartTwo = " + str(ProbabilityPartTwo) + " ProbPartThree = " + str(ProbabilityPartThree) + " ProbPartFour = " + str(ProbabilityPartFour) + " ProbPartFive = " + str(ProbabilityPartFive) + " ProbTotal = " + str(Probability)
	# print "CombinationsPartOne = " + str(CombinVector[0]) + " CombinationsPartTwo = " + str(CombinVector[1]) + " CombinationsPartThree = " + str(CombinVector[2]) + " CombinationsPartFour = " + str(CombinVector[3]) + " CombinationsPartFive = " + str(CombinVector[4])
	return Probability



cdef p_ui_Gi_C_OriginalPre( int PresentDayDerivedAlleleCount, double AlleleFrequency, int AncestralAlleleCounts, int DerivedAlleleCounts, double[:,:] ThisSFS, double ContaminationRate, double[:] BinomVector, int Index, double[:] CombinVector):
	cdef int R_i = AncestralAlleleCounts + DerivedAlleleCounts
	cdef int u_i = DerivedAlleleCounts, j
	cdef int CellToCheck
	# p_kPlus_2 = binom.pmf(2, R_i - u_i, ContaminationRate*AlleleFrequency)
	# p_kPlus_1 = binom.pmf(1, R_i - u_i, ContaminationRate*AlleleFrequency)
	# p_kMinus_1 = binom.pmf(1, u_i, ContaminationRate*( 1 - AlleleFrequency))
	# p_kMinus_2 = binom.pmf(2, u_i, ContaminationRate*( 1 - AlleleFrequency))
	# cdef double p_kPlus_2 = BinomVector[(Index-1)*4]
	# cdef double p_kPlus_1 = BinomVector[(Index-1)*4 + 1]
	# cdef double p_kMinus_1 = BinomVector[(Index-1)*4 + 2]
	# cdef double p_kMinus_2 = BinomVector[(Index-1)*4 + 3]
	# print "Index = " + str(Index)
	# ContaminationRate = ContaminationRate + 1
	cdef int DerCount, AncCount, ThisR_i = AncestralAlleleCounts + DerivedAlleleCounts
	cdef double ProbabilityPartOne = 0, ProbabilityPartTwo = 0, ProbabilityPartThree = 0, ProbabilityPartFour = 0, ProbabilityPartFive = 0, Probability = 0
	# print " Cell to check = " + str(Index) + " Binom = "+ str(BinomVector) + " Pr = " + str(PresentDayDerivedAlleleCount) + "\n"
	for j in range(0, 3):
		CellToCheck = Index + j
		# print "Cell Number = " + str(CellToCheck)
		# Probability = Probability + ThisSFS[j,PresentDayDerivedAlleleCount]
		Probability = Probability + BinomVector[CellToCheck] * ThisSFS[j,PresentDayDerivedAlleleCount]
		# print "Probs = " + str(Probability) + " " + str(BinomVector[CellToCheck]) + " " + str(ThisSFS[j,PresentDayDerivedAlleleCount]) + "\n"
	return Probability

cdef p_ui_Gi_OriginalPre(int DerivedAlleleCounts, int AncestralAlleleCounts, int PresentDayDerivedAlleleCount, double[:,:] ThisSFS, double ContaminationRate, double AlleleFrequency):
	cdef double Probability = 0, TotalSum = 0, FreqWithContamination
	cdef int R_i = DerivedAlleleCounts + AncestralAlleleCounts
	for j in range(0, 3):
		if j == 0:
			FreqWithContamination = ContaminationRate * AlleleFrequency
			Probability = Probability + binom.pmf(DerivedAlleleCounts, R_i, FreqWithContamination) * ThisSFS[j,PresentDayDerivedAlleleCount]
		if j == 1:
			FreqWithContamination = ContaminationRate * AlleleFrequency + (1 - ContaminationRate) * ( 1 / 2)
			Probability = Probability + binom.pmf(DerivedAlleleCounts, R_i, FreqWithContamination) * ThisSFS[j,PresentDayDerivedAlleleCount]
		if j == 2:
			FreqWithContamination = ContaminationRate * AlleleFrequency + (1 - ContaminationRate)
			Probability = Probability + binom.pmf(DerivedAlleleCounts, R_i, FreqWithContamination) * ThisSFS[j,PresentDayDerivedAlleleCount]
	return Probability


cdef ncr(int n, int r):
	if (r > n):
		return 0
	if (r * 2 > n):
		r = n-r
	if (r == 0):
		return 1
	cdef int result, i
	result = n
	for i in range ( 2, r ):
		result *= (n-i+1)
		result /= i
	return result

def p_V_G_given_C_E(int PresentDayDerivedAlleleCount, double AlleleFrequency, int AncestralAlleleCounts, int DerivedAlleleCounts, double EPlus, double EMinus, double EPlusPlus, double EMinusMinus, double EPlusMinus, double ECero, double[:, :] ThisSFS, double ContaminationRate):
	cdef int R_i
	cdef double ProbabilityPartOne, ProbabilityPartTwo, ProbabilityPartThree, ProbabilityPartFour, ProbabilityPartFive, Probability
	R_i = AncestralAlleleCounts + DerivedAlleleCounts
	ProbabilityPartOne = (ECero + EPlusMinus) * p_ui_Gi_C (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts, DerivedAlleleCounts, ThisSFS, ContaminationRate)/ ncr(R_i,DerivedAlleleCounts)
	ProbabilityPartTwo = ( EMinusMinus) * p_ui_Gi_C (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts + 2, DerivedAlleleCounts - 2, ThisSFS, ContaminationRate)/ ncr(R_i,DerivedAlleleCounts - 2)
	ProbabilityPartThree = (EMinus) * p_ui_Gi_C (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts + 1, DerivedAlleleCounts - 1, ThisSFS, ContaminationRate)/ ncr(R_i,DerivedAlleleCounts - 1)
	ProbabilityPartFour = (EPlus) * p_ui_Gi_C (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts - 1, DerivedAlleleCounts + 1, ThisSFS, ContaminationRate)/ ncr(R_i,DerivedAlleleCounts + 1)
	ProbabilityPartFive = (EPlusPlus) * p_ui_Gi_C (PresentDayDerivedAlleleCount, AlleleFrequency, AncestralAlleleCounts - 2, DerivedAlleleCounts + 2, ThisSFS, ContaminationRate)/ ncr(R_i,DerivedAlleleCounts + 2)
	Probability = ProbabilityPartOne + ProbabilityPartTwo + ProbabilityPartThree + ProbabilityPartFour + ProbabilityPartFive
	return Probability

def p_ui_Gi_C( int PresentDayDerivedAlleleCount, double AlleleFrequency, int AncestralAlleleCounts, int DerivedAlleleCounts, double[:, :] ThisSFS, double ContaminationRate):
	cdef int R_i, u_i
	cdef double p_kPlus_2, p_kPlus_1, p_kMinus_1, p_kMinus_2, ProbabilityPartOne, ProbabilityPartTwo, ProbabilityPartThree, ProbabilityPartFour, ProbabilityPartFive, Probability
	R_i = AncestralAlleleCounts + DerivedAlleleCounts
	u_i = DerivedAlleleCounts
	p_kPlus_2 = binom.pmf(2, R_i - u_i, ContaminationRate*AlleleFrequency)
	p_kPlus_1 = binom.pmf(1, R_i - u_i, ContaminationRate*AlleleFrequency)
	p_kMinus_1 = binom.pmf(1, u_i, ContaminationRate*( 1 - AlleleFrequency))
	p_kMinus_2 = binom.pmf(2, u_i, ContaminationRate*( 1 - AlleleFrequency))
	ProbabilityPartOne = p_kPlus_2 * p_ui_Gi(DerivedAlleleCounts + 2, AncestralAlleleCounts - 2,PresentDayDerivedAlleleCount, ThisSFS)
	ProbabilityPartTwo = p_kPlus_1 * p_ui_Gi(DerivedAlleleCounts + 1, AncestralAlleleCounts - 1,PresentDayDerivedAlleleCount, ThisSFS)
	ProbabilityPartThree = (1 - p_kPlus_2 - p_kPlus_1 - p_kMinus_1 - p_kMinus_2) * p_ui_Gi(DerivedAlleleCounts,AncestralAlleleCounts,PresentDayDerivedAlleleCount, ThisSFS)
	ProbabilityPartFour = p_kMinus_1 * p_ui_Gi(DerivedAlleleCounts - 1, AncestralAlleleCounts + 1,PresentDayDerivedAlleleCount, ThisSFS)
	ProbabilityPartFive = p_kMinus_2 * p_ui_Gi(DerivedAlleleCounts - 2, AncestralAlleleCounts + 2,PresentDayDerivedAlleleCount, ThisSFS)
	Probability = ProbabilityPartOne + ProbabilityPartTwo + ProbabilityPartThree + ProbabilityPartFour + ProbabilityPartFive
	return Probability

def p_ui_Gi(int DerivedAlleleCounts, int AncestralAlleleCounts, int PresentDayDerivedAlleleCount, double[:, :] ThisSFS):
	cdef int R_i
	cdef double Probability = 0
	R_i = DerivedAlleleCounts + AncestralAlleleCounts
	for j in range(0, 3):
		if ((PresentDayDerivedAlleleCount == 0) and (j == 0)):
			continue
		if ((PresentDayDerivedAlleleCount == 2) and (j == 2)):
			continue
		if (j==0):
			if (DerivedAlleleCounts == 0):
				Probability = Probability + ThisSFS[PresentDayDerivedAlleleCount,j]
		if (j==1):
			Probability = Probability + ncr(R_i, DerivedAlleleCounts) * (1.0 / (2.0**R_i)) * ThisSFS[PresentDayDerivedAlleleCount,j]
			# print('Here!')
		if (j==2):
			if (AncestralAlleleCounts == 0):
				Probability = Probability + ThisSFS[PresentDayDerivedAlleleCount,j]
	return Probability


