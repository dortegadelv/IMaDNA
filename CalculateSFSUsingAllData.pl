$ProportionsInput = $ARGV[0];
$SFSPrefix = $ARGV[1];
$OutputSFS = $ARGV[2];

open (INPUT,$ProportionsInput) or die "NO!";
@FinalSFS = ();
@PropState = ();
@BranchLengthSingleton = ();
@BranchLengthDoubleton = ();

while (<INPUT>){
chomp;
$Line = $_;
@SplitLine = split (/\s+/, $Line);
push (@PropState,$SplitLine[0]);
push (@BranchLengthSingleton,$SplitLine[1]);
push (@BranchLengthDoubleton,$SplitLine[2]);
}

close (INPUT);

for ($i = 1; $i < 6; $i++){
$SFSFile = $SFSPrefix.".".$i.".txt";
open (SFS,$SFSFile) or die "NO! $SFSFile $SFSPrefix";
$SFSFlag = 0;
@ThisSFS = ();
$RowNumber = 0;
while (<SFS>){
chomp;
$Line = $_;
if ($Line =~ "Expected AFS"){
$SFSFlag = 0;
}

if ($SFSFlag == 1){

@LineElements = split (/\s+/, $Line);
for ($j = 0; $j < scalar(@LineElements); $j++){
$ThisSFS[$RowNumber][$j] = $LineElements[$j];
}
$ColumnNumber = scalar(@LineElements);
$RowNumber++;
}

if ($Line =~ "Matrix sum"){
# print "Here!\n";
$SFSFlag = 1;
}
}
close (SFS);

### Cases

## Case 1
if ($i == 1){
@ReformatBranchLengths = ();
# print "Branch lengths before formatting\n";
for ($k = 0; $k < $RowNumber; $k++){
for ($j = 0; $j < scalar(@LineElements); $j++){
# print "$ThisSFS[$k][$j]\t";
}
# print "\n";
}

# print "Branch lengths\n";
$BranchLengthSingletonSum = 0;
for ($k = 0; $k < $RowNumber; $k++){
for ($j = 0; $j < scalar(@LineElements); $j++){
if ( ( $k == 0 ) && ( $j == 0 )){
$ReformatBranchLengths[$j][$k] = 0;
} else {
$ReformatBranchLengths[$j][$k] = $ThisSFS[$k][$j] - $PastBL;
}
if ( $k == 1 ){
$BranchLengthSingletonSum = $BranchLengthSingletonSum + $ThisSFS[$k][$j];
}
$PastBL = $ThisSFS[$k][$j];
# print "$ReformatBranchLengths[$k][$j]\t";
}
# print "\n";
}
# print "End print branch lengths\n";
$ReformatBranchLengths[0][1] = $ReformatBranchLengths[0][1] + $BranchLengthSingleton[0];

$TotalSum = 0;
for ($k = 0; $k < $RowNumber; $k++){
for ($j = 0; $j < scalar(@LineElements); $j++){
$TotalSum = $TotalSum + $ReformatBranchLengths[$k][$j];
}
# print "\n"
}

for ($k = 0; $k < $RowNumber; $k++){
for ($j = 0; $j < scalar(@LineElements); $j++){
$ReformatBranchLengths[$k][$j] = $ReformatBranchLengths[$k][$j] ;
}
# print "\n";
}

# print "Exit\n";

# print "After case 1 SFS\n";
# print "Proportion 1 = $PropState[$i]\n";
for ($k = 0; $k < $RowNumber; $k++){
for ($j = 0; $j < scalar(@LineElements); $j++){
# print "$ReformatBranchLengths[$k][$j]\t";
$FinalSFS[$k][$j] = $ReformatBranchLengths[$k][$j] * $PropState[$i-1];
# print "$FinalSFS[$k][$j]\t";
}
# print "\n";
}
}

## Case 2
if ($i == 2){
@ReformatBranchLengths = ();
$ReformatBranchLengths[0][1] = $BranchLengthSingleton[2] + ($ThisSFS[1][0] - $ThisSFS[0][3]) + ($ThisSFS[0][1] - $ThisSFS[0][0])/3 ;
$ReformatBranchLengths[1][1] = 2*($ThisSFS[1][1] - $ThisSFS[1][0])/3 + 2*($ThisSFS[0][2] - $ThisSFS[0][1])/3;
$ReformatBranchLengths[2][1] = ($ThisSFS[0][3] - $ThisSFS[0][2]) + ($ThisSFS[1][2] - $ThisSFS[1][1])/3 ;
$ReformatBranchLengths[0][0] = 0;
$ReformatBranchLengths[1][0] = 2*($ThisSFS[0][1] - $ThisSFS[0][0])/3;
$ReformatBranchLengths[2][0] = ($ThisSFS[0][2] - $ThisSFS[0][1])/3;
$ReformatBranchLengths[0][2] = ($ThisSFS[1][1] - $ThisSFS[1][0])/3;
$ReformatBranchLengths[1][2] = 2* ($ThisSFS[1][2] - $ThisSFS[1][1])/3;
$ReformatBranchLengths[2][2] = 0 ;

$TotalSum = 0;
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$TotalSum = $TotalSum + $ReformatBranchLengths[$k][$j];
}
# print "\n"
}

for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$ReformatBranchLengths[$k][$j] = $ReformatBranchLengths[$k][$j] ;
}
# print "\n";
}

# print "After case 2 SFS\n";
# print "Proportion 3 = $PropState[$i]\n";
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
# print "$ReformatBranchLengths[$k][$j]\t";
$FinalSFS[$k][$j] = $FinalSFS[$k][$j] + $ReformatBranchLengths[$k][$j] * $PropState[$i-1];
# print "$FinalSFS[$k][$j]\t";
}
# print "\n";
}

}

## Case 3
if ($i == 3){
@ReformatBranchLengths = ();
$ReformatBranchLengths[0][1] = $BranchLengthSingleton[2] + ($ThisSFS[0][1] - $ThisSFS[0][0] + $ThisSFS[1][1] - $ThisSFS[1][0])/2;
$ReformatBranchLengths[1][1] = 4*($ThisSFS[0][2] - $ThisSFS[0][1] + $ThisSFS[1][2] - $ThisSFS[1][1])/6;
$ReformatBranchLengths[2][1] = ($ThisSFS[0][3] - $ThisSFS[0][2] + $ThisSFS[1][3] - $ThisSFS[1][2])/2;
$ReformatBranchLengths[0][0] = 0;
$ReformatBranchLengths[1][0] = ($ThisSFS[0][1] - $ThisSFS[0][0] + $ThisSFS[1][1] - $ThisSFS[1][0])/2;
$ReformatBranchLengths[2][0] = ($ThisSFS[0][2] - $ThisSFS[0][1] + $ThisSFS[1][2] - $ThisSFS[1][1])/6;
$ReformatBranchLengths[0][2] = ($ThisSFS[0][2] - $ThisSFS[0][1] + $ThisSFS[1][2] - $ThisSFS[1][1])/6;
$ReformatBranchLengths[1][2] = ($ThisSFS[0][3] - $ThisSFS[0][2] + $ThisSFS[1][3] - $ThisSFS[1][2])/2;
$ReformatBranchLengths[2][2] = 0 ;

$TotalSum = 0;
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$TotalSum = $TotalSum + $ReformatBranchLengths[$k][$j];
}
# print "\n"
}

for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$ReformatBranchLengths[$k][$j] = $ReformatBranchLengths[$k][$j] ;
}
# print "\n";
}

# print "After case 3 SFS\n";
# print "Proportion 3 = $PropState[$i]\n";
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
# print "$ReformatBranchLengths[$k][$j]\t";
$FinalSFS[$k][$j] = $FinalSFS[$k][$j] + $ReformatBranchLengths[$k][$j] * $PropState[$i-1];
# print "$FinalSFS[$k][$j]\t";
}
# print "\n";
}


}

## Case 4
if ($i == 4){

@ReformatBranchLengths = ();
for ($k = 0; $k < $RowNumber; $k++){
for ($j = 0; $j < scalar(@LineElements); $j++){
# print "$ThisSFS[$k][$j]\t";
}
# print "\n";
}
$ReformatBranchLengths[0][1] = $BranchLengthSingleton[3];
$ReformatBranchLengths[1][1] = 0;
$ReformatBranchLengths[2][1] = 0;
$ReformatBranchLengths[0][0] = $ThisSFS[0][0];
$ReformatBranchLengths[1][0] = ($ThisSFS[0][1] - $ThisSFS[0][0]);
$ReformatBranchLengths[2][0] = $ThisSFS[0][2] - $ThisSFS[0][1];
$ReformatBranchLengths[0][2] = $ThisSFS[1][0] - $ThisSFS[0][2];
$ReformatBranchLengths[1][2] = $ThisSFS[1][1] - $ThisSFS[1][0];
$ReformatBranchLengths[2][2] = $ThisSFS[1][2] - $ThisSFS[1][1] ;

$DoubletonSum = $ReformatBranchLengths[0][2] + $ReformatBranchLengths[1][2] + $ReformatBranchLengths[2][2];
$ReformatBranchLengths[0][2] = $ReformatBranchLengths[0][2] + $BranchLengthDoubleton[3];
# $ReformatBranchLengths[1][2] = $ReformatBranchLengths[1][2] + $BranchLengthDoubleton[3] * $ReformatBranchLengths[1][2] / $DoubletonSum;
# $ReformatBranchLengths[2][2] = $ReformatBranchLengths[2][2] + $BranchLengthDoubleton[3] * $ReformatBranchLengths[2][2] / $DoubletonSum;

# print "Refformatted data\n";
$TotalSum = 0;
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$TotalSum = $TotalSum + $ReformatBranchLengths[$k][$j];
}
# print "\n";
}

for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$ReformatBranchLengths[$k][$j] = $ReformatBranchLengths[$k][$j];
$FinalSFS[$k][$j] = $FinalSFS[$k][$j] + $ReformatBranchLengths[$k][$j] * $PropState[$i-1];
}
# print "\n";
}
}

## Case 5
if ($i == 5){
@ReformatBranchLengths = ();

$ReformatBranchLengths[0][1] = $BranchLengthSingleton[4];
$ReformatBranchLengths[1][1] = 0;
$ReformatBranchLengths[2][1] = 0;
$ReformatBranchLengths[0][0] = 0;
$ReformatBranchLengths[1][0] = ($ThisSFS[0][1] - $ThisSFS[0][0] + $ThisSFS[1][1] - $ThisSFS[1][0])*2/3;
$ReformatBranchLengths[2][0] = ($ThisSFS[0][2] - $ThisSFS[0][1] + $ThisSFS[1][2] - $ThisSFS[1][1])/3;
$ReformatBranchLengths[0][2] = $BranchLengthDoubleton[3] + ($ThisSFS[0][1] - $ThisSFS[0][0] + $ThisSFS[1][1] - $ThisSFS[1][0])/3;
$ReformatBranchLengths[1][2] = ($ThisSFS[0][2] - $ThisSFS[0][1] + $ThisSFS[1][2] - $ThisSFS[1][1])*2/3;
$ReformatBranchLengths[2][2] = 0;

$TotalSum = 0;
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$TotalSum = $TotalSum + $ReformatBranchLengths[$k][$j];
}
# print "\n";
}

for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$ReformatBranchLengths[$k][$j] = $ReformatBranchLengths[$k][$j];
$FinalSFS[$k][$j] = $FinalSFS[$k][$j] + $ReformatBranchLengths[$k][$j] * $PropState[$i-1];
}
# print "\n";
}


}

}

# print "Final SFS\n\n";
$TotalSum = 0;
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
$TotalSum = $TotalSum + $FinalSFS[$k][$j];
# $FinalSFS[$k][$j] = $FinalSFS[$k][$j] / $TotalSum;
}
}

# $TotalSum = 0;
for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
# $TotalSum = $TotalSum + $FinalSFS[$k][$j];
$FinalSFS[$k][$j] = $FinalSFS[$k][$j] / $TotalSum;
}
}


open (OUTPUT,">$OutputSFS") or die "NO!";

for ($k = 0; $k < 3; $k++){
for ($j = 0; $j < 3; $j++){
print OUTPUT "$FinalSFS[$k][$j]\t";
}
print OUTPUT "\n";
}

close (OUTPUT);

