#### If we sample two chromosomes from the present, there are 5 possible ways the chromosomes
#### can be present at the time that we sample two chromosomes from an ancestral individual

NeOne=$1
NeTwo=$2
NeAncestral=$3
MutationRate=1
RecRate=0.00000001
MigRateOneToTwo=$5
MigRateTwoToOne=$6
DivTime=$4
AncTime=$7

DivTime=$( echo "$DivTime - $AncTime" | bc -l  )
ThetaOne=$( echo "4 * $NeOne" | bc -l )
# echo $ThetaOne
ThetaTwo=$( echo "4 * $NeTwo" | bc -l )
# echo $ThetaTwo
ThetaAncestral=$( echo "4 * $NeAncestral" | bc -l )
# echo $ThetaAncestral
FractionN2=$( echo "$NeTwo / $NeOne" | bc -l )
FractionNa=$( echo "$NeAncestral / $NeOne" | bc -l )
# echo $FractionN2
# echo $FractionNa
MsDivTime=$( echo "$DivTime / ( 4 * $NeOne)" | bc -l )
MigRateOneToTwo=$( echo $MigRateOneToTwo | awk '{ print sprintf("%.18f", $1); }')
# MigRateOneToTwo=$( echo "4 * $NeOne * $MigRateOneToTwo" | bc -l )
mOne=$( echo "$MigRateOneToTwo / ( 4 * $NeOne )" | bc -l )
MigRateTwoToOne=$( echo $MigRateTwoToOne | awk '{ print sprintf("%.18f", $1); }')
# MigRateTwoToOne=$( echo "4 * $NeOne * $MigRateTwoToOne" | bc -l )
mTwo=$( echo "$MigRateTwoToOne / ( 4 * $NeOne )" | bc -l )
MsAncTime=$( echo "$AncTime / $NeOne" | bc -l )
OutputProportions=$9
PrefixIMCLAM=$8
OutputFinalSFS=${10}
DivisionsInCoalescentCalculations=${11}
rm $OutputProportions

# echo "PrintParams"
# echo $NeOne
# echo $NeTwo
# echo $NeAncestral
# echo $DivTime
# echo $MigRateOneToTwo
# echo $MigRateTwoToOne
# echo $AncTime
# echo $PrefixIMCLAM
# echo $OutputProportions
# echo $OutputFinalSFS

#### 1) The two present day chromosomes are in population 1

OutputSFS=$PrefixIMCLAM."1.txt"
/opt/IMaDNA/im_clam -s /opt/IMaDNA/testSS_22 -m /opt/IMaDNA/testSS_22_mats -exp -x $FractionN2,$FractionNa,$MigRateOneToTwo,$MigRateTwoToOne,$MsDivTime > $OutputSFS
# Rscript --vanilla CoalTimePresentBeforeAncestral.R $mOne $mTwo $ThetaOne $ThetaTwo $AncTime 1000 1 $MutationRate $OutputProportions
echo "$PETSC_DIR/arch-linux2-c-debug/bin/mpiexec -n 1 ../../../im_clam-master/im_clam -s ../../../im_clam-master/stateSpaceFiles/testSS_22 -m ../../../im_clam-master/stateSpaceFiles/testSS_22_mats -exp -x $FractionN2,$FractionNa,$MigRateOneToTwo,$MigRateTwoToOne,$MsDivTime"
#### 2) One present day chromosome is in population 1, the other one is in population 2

OutputSFS=$PrefixIMCLAM."2.txt"
/opt/IMaDNA/im_clam -s /opt/IMaDNA/ss_1_3 -m /opt/IMaDNA/ss_1_3_mats -exp -x $FractionN2,$FractionNa,$MigRateOneToTwo,$MigRateTwoToOne,$MsDivTime > $OutputSFS
# Rscript --vanilla CoalTimePresentBeforeAncestral.R $mOne $mTwo $ThetaOne $ThetaTwo $AncTime 1000 2 $MutationRate $OutputProportions

#### 3) The two present day chromosomes are in population 2

OutputSFS=$PrefixIMCLAM."3.txt"
/opt/IMaDNA/im_clam -s /opt/IMaDNA/ss_1_4 -m /opt/IMaDNA/ss_1_4_mats -exp -x $FractionN2,$FractionNa,$MigRateOneToTwo,$MigRateTwoToOne,$MsDivTime > $OutputSFS
# Rscript --vanilla CoalTimePresentBeforeAncestral.R $mOne $mTwo $ThetaOne $ThetaTwo $AncTime 1000 3 $MutationRate $OutputProportions

#### 4) The two present day chromosomes have coalesced, and that lineage is now present in population 1

OutputSFS=$PrefixIMCLAM."4.txt"
/opt/IMaDNA/im_clam -s /opt/IMaDNA/ss_1_2 -m /opt/IMaDNA/ss_1_2_mats -exp -x $FractionN2,$FractionNa,$MigRateOneToTwo,$MigRateTwoToOne,$MsDivTime > $OutputSFS
# Rscript --vanilla CoalTimePresentBeforeAncestral.R $mOne $mTwo $ThetaOne $ThetaTwo $AncTime 1000 4 $MutationRate $OutputProportions

#### 5) The two present day chromosomes have coalesced, and that lineage is now present in population 2

OutputSFS=$PrefixIMCLAM."5.txt"
/opt/IMaDNA/im_clam -s /opt/IMaDNA/ss_1_3 -m /opt/IMaDNA/ss_1_3_mats -exp -x $FractionN2,$FractionNa,$MigRateOneToTwo,$MigRateTwoToOne,$MsDivTime > $OutputSFS
# Rscript --vanilla CoalTimePresentBeforeAncestral.R $mOne $mTwo $ThetaOne $ThetaTwo $AncTime 1000 5 $MutationRate $OutputProportions
Rscript --vanilla /opt/IMaDNA/CoalTimePresentBeforeAncestral.R $mOne $mTwo $ThetaOne $ThetaTwo $AncTime $DivisionsInCoalescentCalculations 6 $MutationRate $OutputProportions


# echo $OutputProportions
# echo $PrefixIMCLAM
# echo $OutputFinalSFS

perl /opt/IMaDNA/CalculateSFSUsingAllData.pl $OutputProportions $PrefixIMCLAM $OutputFinalSFS


