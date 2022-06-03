import numpy as np
# Initialize matrix main diagonal and diagonal below main diagonal and with zeros
#Input : Matrix size
#Output : Matrix filled with nan while main diagonal and diagonal below main diagonal filled with zeros

def InitializeMatrix(N):
    Matrix = np.empty((N, N))
    Matrix.fill(np.nan)
    for i in range(len(Matrix)):
        for j in range(len(Matrix[i])):
            if i == j:
                Matrix[i][j] = 0
                if i == 1:
                    continue
                Matrix[i - 1][j - 2] = 0

    return Matrix

# Get maximum value of K from Bifurcation (merging two substructures)
# This should be called in the fill matrix function
# This should return the value of k that gives the max value and the index of k

def bifurication(i, j, matrix):
    valuesOfK = []
    maxK = 0
    index = 0
    for k in range(i, j):
        if k > i and k < j:
            valuesOfK.append(k)
    for k in valuesOfK:
        maxK = max(matrix[i, k] + matrix[k + 1, j], maxK)
        if maxK == matrix[i, k] + matrix[k + 1, j]:
            index = k
    return maxK, index
   #return maxK, maxKIndex


# Check if matrix[i] and matrix[j] are paired or not
# Input : takes bases from function fill matrix
# Output : 1 or 0
def diagonal_check(bases):

    if (bases[0] == 'A' and bases[1]=='U' or  bases[0] == 'U' and  bases[1] =='A')or(bases[0] == 'G' and  bases[1]=='C' or  bases[0] == 'C' and  bases[1] =='G')or(bases[0] == 'G' and  bases[1]=='U' or  bases[0] == 'U' and  bases[1] =='G'):
        return 1
    else:
        return 0
    #return 1 if bases are paired else return 0


# Mapping each character in the input sequence to a unique identifier
# to access them when filling the matrix
# Input : RNA sequence
# Output : dictionary containing the sequence mapped to a unique key
# example seq : GGGAAAUCC , dictionary : {0:G,1:G,2:G,3:A,4:A,5:A,6:U,7:C,8:C}
def rnaDictionary(sequence):
    seqList = list(sequence)
    dict = {}
    for i in range(len(sequence)):
        dict[i] = seqList[i]
    return dict
    #return rnaDict

#Should fill the matrix
# Input  : matrix , rnadictionary, size of matrix
# no return as the original matrix is modified
# Mustafa and Ibrahim
def FillMatrix(matrix, rnaDict, N) :
    for m in range(1, N):
        for i in range(0, N - m):
            j = i + m
            bases = (rnaDict[i], rnaDict[j])
            down = matrix[i + 1][j]
            left = matrix[i][j - 1]
            diagonal = matrix[i + 1][j - 1] + diagonal_check(bases)
            k, index = bifurication(i, j, matrix)
            matrix[i][j] = max(down, left, diagonal, k)


#Should iterate over the original matrix and return the final list
def traceback(row,column,matrix,outputArr,rnaDict):
    if(row>column or row==column):
        return outputArr

    while True:
        currVal = matrix[row][column]
        if (currVal == 0):
            break
        down = matrix[row + 1][column]
        left = matrix[row][column - 1]
        diagonal = matrix[row + 1][column - 1]
        bases = (rnaDict[row], rnaDict[column])
        if (currVal == down):
            row = row + 1
            continue

        elif (currVal == left):
            column = column - 1
            continue

        elif (currVal == diagonal + diagonal_check(bases)):
            if (row < column):
                outputArr[row] = "("
                outputArr[column] = ")"
            else:
                outputArr[column] = "("
                outputArr[row] = ")"
            row = row + 1
            column = column - 1

            continue

        else:
            ThemaxK,TheK= bifurication(row, column, matrix)
            traceback(row,TheK,matrix,outputArr,rnaDict)
            traceback(TheK,column,matrix,outputArr,rnaDict)
            break
    return outputArr
    #return finalList

#Main
if __name__ == "__main__":
    rna = input("please input the RNA Sequence: ")
    rna = rna.upper()
    RnaMatrix = InitializeMatrix(len(rna))
    rnaDict = rnaDictionary(rna)
    FillMatrix(RnaMatrix, rnaDict, len(rna))
    TheoutputArray = ['.'] * len(rna)
    column = len(rna) - 1
    row = 0
    TheoutputArray = traceback(row, column, RnaMatrix, TheoutputArray, rnaDict)
    print(''.join(TheoutputArray))