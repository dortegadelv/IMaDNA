
library(expm)

args = commandArgs(trailingOnly=TRUE)
# print (args[1])
# print (args[2])
# print (args[3])
# print (args[4])
# print (args[5])
# print (args[6])
# print (args[7])
Matrix <- matrix(,nrow=5,ncol=5)

Value = as.numeric(args[1]) + as.numeric(args[2])
# print (Value)

m1 <- as.numeric(args[1])
m2 <- as.numeric(args[2])
Theta1 <- as.numeric(args[3])
Theta2 <- as.numeric(args[4])
ThetaA <- 0.002
DivTime <- 0.003

#################### Set mu ###########################

mu <- as.numeric(args[8])
N_1 <- Theta1 / (4*mu)
N_2 <- Theta2 / (4*mu)
N_A <- ThetaA / (4*mu)
mig1 <- m1 * (mu)
mig2 <- m2 * (mu)

FourN1mig1 <- 4*N_1* mig1
FourN2mig2 <- 4*N_2* mig2
DivergenceTime <- as.numeric(args[5])
DivTimeMs <- DivergenceTime * mu

ActualTimeMs <- DivergenceTime / (4 * N_1)
# DivTimeMs <- as.numeric(args[5])
DivisionNumber <- as.numeric(args[6])
TimeDivisions <- DivTimeMs / as.numeric ( args[6] )

OutputProportionsFile <- args[9]
OutputVector <- c()
#######################################################


Matrix[1,1] <- -2*m1 - 2/Theta1
Matrix[1,2] <- 2*m1
Matrix[1,3] <- 0
Matrix[1,4] <- 2/Theta1
Matrix[1,5] <- 0

Matrix[2,1] <- m2
Matrix[2,2] <- -m2 - m1
Matrix[2,3] <- m1
Matrix[2,4] <- 0
Matrix[2,5] <- 0

Matrix[3,1] <- 0
Matrix[3,2] <- 2*m2
Matrix[3,3] <- -2*m2 - 2/Theta2
Matrix[3,4] <- 0
Matrix[3,5] <- 2/Theta2

Matrix[4,1] <- 0
Matrix[4,2] <- 0
Matrix[4,3] <- 0
Matrix[4,4] <- -m1
Matrix[4,5] <- m1

Matrix[5,1] <- 0
Matrix[5,2] <- 0
Matrix[5,3] <- 0
Matrix[5,4] <- m2
Matrix[5,5] <- -m2

AllFT_11 <- c()
AllFT_12 <- c()
AllFT_22 <- c()

TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
if ((args[7] == 1) | (args[7] == 6)){
TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
# print("Ending Probability")
# print(MatrixExponentiation[1,1])
# print ("Prob of singleton mutation")
# print (NoMutation)
ActTime <- DivTimeMs / (mu)
ActTime <- ActTime / (4 * N_1)
# print ("Length of singleton branch")
BranchLength = 2 * ActTime
# print (BranchLength)
OutputVector <- c(MatrixExponentiation[1,1],BranchLength,0)
write.table(t(OutputVector), file = OutputProportionsFile, sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

if (args[7] == 2 | (args[7] == 6)){
TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
# print("Ending Probability")
# print(MatrixExponentiation[1,2])
# print ("Prob of singleton mutation")
# print (NoMutation)
ActTime <- DivTimeMs / (mu)
ActTime <- ActTime / (4 * N_1)
# print ("Length of singleton branch")
BranchLength = 2 * ActTime
# print (BranchLength)
OutputVector <- c(MatrixExponentiation[1,2],BranchLength,0)
write.table(t(OutputVector), file = OutputProportionsFile, sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

if (args[7] == 3 | (args[7] == 6)){
TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
# print("Ending Probability")
# print(MatrixExponentiation[1,3])
# print ("Prob of singleton mutation")
# print (NoMutation)
ActTime <- DivTimeMs / (mu)
ActTime <- ActTime / (4 * N_1)
# print ("Length of singleton branch")
BranchLength = 2 * ActTime
# print (BranchLength)
OutputVector <- c(MatrixExponentiation[1,3],BranchLength,0)
write.table(t(OutputVector), file = OutputProportionsFile, sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
Sum <- MatrixExponentiation[1,1] + MatrixExponentiation[1,2] + MatrixExponentiation[1,3]
# print ("The sum of the three components = ")
# print (Sum)
Sum <- MatrixExponentiation[1,4] + MatrixExponentiation[1,5]
# print ("And the other two components = ")
# print (MatrixExponentiation[1,4])
# print (MatrixExponentiation[1,5])

if (args[7] == 4 | (args[7] == 6)){
TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
# print("Ending Probability")
# print(MatrixExponentiation[1,4])
OutputVector <- c(MatrixExponentiation[1,4])
ProbNoMutation <- 0
ProbNoDoubleton <- 0
ProbSum <- 0
for (i in 1:DivisionNumber){
Time <- i * TimeDivisions - TimeDivisions / 2
TimeMatrix <- Matrix * Time
MatrixExponentiation <- expm(TimeMatrix)
OtherTime <- (DivisionNumber - i) * TimeDivisions + TimeDivisions / 2
TimeMatrix <- Matrix * OtherTime
MatrixExponentiationTwo <- expm(TimeMatrix)
ProbSum <- ProbSum + ( MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiationTwo[4,4] + MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiationTwo[5,4] )
}

ProbOne <- 0
for (i in 1:DivisionNumber){
Time <- i * TimeDivisions - TimeDivisions / 2
TimeMatrix <- Matrix * Time
MatrixExponentiation <- expm(TimeMatrix)
OtherTime <- (DivisionNumber - i) * TimeDivisions + TimeDivisions / 2
TimeMatrix <- Matrix * OtherTime
MatrixExponentiationTwo <- expm(TimeMatrix)
F1_11 <- MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiation[4,4] * TimeDivisions
F2_11 <- MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiation[5,4] * TimeDivisions
FT_11 <- F1_11 + F2_11
CurTime <- (i - 0.5) / DivisionNumber * DivTimeMs
CurTime <- CurTime / (mu)
CurTime <- CurTime / (4 * N_1)
ProbNoMutation <- ProbNoMutation + ( MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiationTwo[4,4] + MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiationTwo[5,4] ) / ProbSum * 2*CurTime
ProbNoDoubleton <- ProbNoDoubleton + ( MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiationTwo[4,4] + MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiationTwo[5,4] ) / ProbSum * (ActualTimeMs - CurTime)
AllFT_11 <- c(AllFT_11,FT_11)
}
# print ("Current time")
# print (CurTime)
# print ("Length of singleton branch")
# print (ProbNoMutation)
# print ("Length of doubleton branch")
# print (ProbNoDoubleton)
# print ("Cur Time")
# print (CurTime)
# print ("Ms time")
# print (DivTimeMs)
OutputVector <- c(OutputVector,ProbNoMutation,ProbNoDoubleton)
write.table(t(OutputVector), file = OutputProportionsFile, sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

if (args[7] == 5 | (args[7] == 6)){
TimeMatrix <- Matrix * DivTimeMs
MatrixExponentiation <- expm(TimeMatrix)
# print("Ending Probability")
# print(MatrixExponentiation[1,5])
OutputVector <- c(MatrixExponentiation[1,5])
ProbNoMutation <- 0
ProbNoDoubleton <- 0
ProbSum <- 0
for (i in 1:DivisionNumber){
Time <- i * TimeDivisions - TimeDivisions / 2
TimeMatrix <- Matrix * Time
MatrixExponentiation <- expm(TimeMatrix)
OtherTime <- (DivisionNumber - i) * TimeDivisions + TimeDivisions / 2
TimeMatrix <- Matrix * OtherTime
MatrixExponentiationTwo <- expm(TimeMatrix)
ProbSum <- ProbSum + ( MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiationTwo[4,5] + MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiationTwo[5,5] )
}

for (i in 1:DivisionNumber){
Time <- i * TimeDivisions - TimeDivisions / 2
TimeMatrix <- Matrix * Time
MatrixExponentiation <- expm(TimeMatrix)
OtherTime <- (DivisionNumber - i) * TimeDivisions + TimeDivisions / 2
TimeMatrix <- Matrix * OtherTime
MatrixExponentiationTwo <- expm(TimeMatrix)
F1_11 <- MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiation[4,5] * TimeDivisions
F2_11 <- MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiation[5,5] * TimeDivisions
FT_11 <- F1_11 + F2_11
CurTime <- (i - 0.5) / DivisionNumber * DivTimeMs
CurTime <- CurTime / (mu)
CurTime <- CurTime / (4 * N_1)
ProbNoMutation <- ProbNoMutation + ( MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiationTwo[4,5] + MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiationTwo[5,5] ) / ProbSum * 2 *CurTime
ProbNoDoubleton <- ProbNoDoubleton + ( MatrixExponentiation[1,1]*(2/Theta1)*MatrixExponentiationTwo[4,5] + MatrixExponentiation[1,3]*(2/Theta2)*MatrixExponentiationTwo[5,5] ) / ProbSum * (ActualTimeMs - CurTime)
AllFT_11 <- c(AllFT_11,FT_11)
}
# print ("Length of singleton branch")
# print (ProbNoMutation)
# print ("Length of doubleton branch")
# print (ProbNoDoubleton)
OutputVector <- c(OutputVector,ProbNoMutation,ProbNoDoubleton)
write.table(t(OutputVector), file = OutputProportionsFile, sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}


