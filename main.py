##############################################
# ############### 1st method ############### #
def LUdecomposition(mat, b):
    """
    :param mat: the coefficients matrix
    :param b:  the solution vector
    :return: none
    """
    inversL, L = toUpperTriangularMat(mat)  # calculate the L matrix and the inverse L
    print("L = ")
    printMatrix(L)  # print the L matrix
    U = multMat(inversL, mat)  # calculate the U matrix
    print("U = ")
    printMatrix(U)  # print the U matrix
    inversU = FromUpperToInvers(U, createIdentityMatrix(len(U)))  # calculate thr inverse of U
    x = multMat(inversL, b)  # finding the result vector
    x = multMat(inversU, x)  # finding the result vector
    print("X = ")
    printMatrix(x)  # print the X matrix


def FromUpperToInvers(A, inverseMat):
    """
    :param A: upper matrix
    :param inverseMat: the matrix that will become the inverse
    :return: Inverse matrix
    """
    elemntarMat = createIdentityMatrix(len(A))  # identity matrix
    for i in range(len(A) - 1, -1, -1):  # run over the columns
        for j in range(i):  # run over the lines above the pivot
            elemntarMat[j][i] = -(A[j][i] / A[i][i])
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            elemntarMat[j][i] = 0
        if A[i][i] != 1:  # convert the pivots to one
            elemntarMat[i][i] = 1 / A[i][i]
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            elemntarMat[i][i] = 1
    return inverseMat


def toUpperTriangularMat(A):
    """
    :param A: the matrix we want to turn into a upper triangle matrics
    :return: the multiplication of the elementary matrics that create U
    """
    L = createIdentityMatrix(len(A))  # create indetity matrix
    InL = createIdentityMatrix(len(A))  # create indetity matrix
    for i in range(len(A)):  # run over the lines
        if A[i][i] == 0:  # if the pivot is 0
            for j in range(i + 1, len(A)):  # run over the columns
                if A[j][i] != 0:  # if the element under the pivot is not 0
                    L = multMat((L, linesExchange(A, i, j)))  # make lines exchange and multiply
                    InL = multMat((linesExchange(A, i, j)), InL)  # make lines exchange and multiply
                    A = multMat((linesExchange(A, i, j)), A)
                    break
        if A[i][i] != 0:  # check if B is regular
            for j in range(i + 1, len(A)):  # run over the columns
                identity = createIdentityMatrix(len(A))
                identity[j][i] = -(A[j][i] / A[i][i])  # elementary matrix
                InL = multMat(identity, InL)  # L^(-1)
                A = multMat(identity, A)
                identity[j][i] *= -1  # changing the element in order to find L
                L = multMat(L, identity)
    return InL, L


def linesExchange(A, line1, line2):
    """
    :param A: A matrix
    :param line1: A line
    :param line2: The line we want to exchange with
    :return: elementry matrix
    """
    idendityMax = createIdentityMatrix(len(A))  # create identity matrix

    # exchange the members in line1
    temp = idendityMax[line1][line1]
    idendityMax[line1][line1] = idendityMax[line2][line1]
    idendityMax[line2][line1] = temp

    # exchange the members in line2
    temp = idendityMax[line2][line2]
    idendityMax[line2][line2] = idendityMax[line1][line2]
    idendityMax[line1][line2] = temp
    return idendityMax


def createIdentityMatrix(size):
    """
    :param size: the size of the square matrix
    :return:
    """
    identityMat = newMat(size, size)  # create a zero matrix in the required size
    for index in range(size):  # go over the main diagonal
        identityMat[index][index] = 1  # change the elements in the main diagonal to 1
    return identityMat


def multMat(A, B):
    """
    :param A: a matrix in sise n*m
    :param B: a mtrix in size m*k
    :return: A*B  (in size n*k)
    """
    if len(A[1]) == len(B):  # check if A and B have the same number of rows and columns
        C = newMat(len(A), len(B[0]))  # the matrix the function returns
        for i in range(len(C)):
            for j in range(len(C[1])):
                for k in range(len(B)):
                    C[i][j] += A[i][k] * B[k][j]
        return C
    else:
        return None  # the multiplication  is impossible


# ############### 1st method ############### #

##############################################

# ############### 2nd method ############### #
def iterativGuassSeidel(A, b, epsilon, flagD):
    """
    :param A: a matrix
    :param b: the result vector
    :param epsilon: the mistake
    :param flagD: tell us if the system have dominant diagonal
    :return: None
    """
    print("** Iterative Guass Seidel **")
    flagC = False  # flagC = false if the linear equations does not converge
    x = newMat(len(b), 1)  # create the result vector
    print("The results are:\nThe first guess is: ")
    printMatrix(x)
    for _ in range(99):  # max number of iterations is 99
        oldX1 = x[0][0]  # copy the old x value of the current iteration
        for i in range(len(x)):  # go over the all variables
            if A[i][i] == 0:  # preventing division by zero
                return None
            temp = b[i][0] / A[i][i]
            for j in range(len(x)):  # calculate the value of the variable in the new iteration
                if i != j:
                    temp -= (A[i][j] * x[j][0]) / A[i][i]
            x[i][0] = temp  # update the calculated value
        print("The result of the " + str(_ + 1) + " iteration is: ")
        printMatrix(x)
        if abs(oldX1 - x[0][0]) < epsilon:  # check stop condition
            flagC = True
            break
    if flagC is True:
        print("***********")  # print final results
        if flagD is False:
            print("Although there is no dominant diagonal, the results are: \n")
        print("The final result is: x = ")
        printMatrix(x)
        print("The number of the iteration is : " + str(_ + 1))
    else:
        print("The system of linear equations does not converge :(")


# ############### 2nd method ############### #


def printMatrix(a):
    """
    :param a: a matrix to print
    :return: prints in matrix format
    """
    print("-----------------------------")
    for i in range(len(a)):
        if i is len(a) - 1:
            print(" " + str(a[i]) + "]")
        elif i is 0:
            print("[" + str(a[i]))
        else:
            print(" " + str(a[i]))
    print("-----------------------------")


def newMat(numRow, numCol):
    """
    :param numRow: the number of rows in the mat
    :param numCol: the number of columns in the mat
    :return: a zero matrix in the required size
    """
    mat = []  # the zero matrix the function returns
    for i in range(numRow):
        mat.append([])  # create a new row
        for j in range(numCol):
            mat[i].append(0)  # fill the row with
    return mat


def dominantDiagonal(A):
    """
    :param A: a list - matrix
    :return: true if the matrix have dominant diagonal else return false
    """
    for i in range(len(A)):
        lineSum = 0  # sum of every line except to the element in the diagonal
        for j in range(len(A)):
            if i != j:
                lineSum += A[i][j]
        if A[i][i] <= lineSum:
            # print("There is no dominant diagonal ...")
            return False
    print("There is a dominant diagonal :)")
    return True


# dominant diagonal part

def copyMat(A):
    """
    :param A: a matrix
    :return: a copy of the matrix
    """
    B = newMat(len(A), len(A[0]))  # create a zero matrix of the same size as A
    for i in range(len(A)):  # copy A
        for j in range(len(A[0])):
            B[i][j] = A[i][j]
    return B


def createDominantDiagonal(A, b=None):
    """
    :param A: a coefficients matrix
    :param b: the column vector of constant terms.
    :return: matrix A with dominant diagonal
    """
    max = 0  # the max value in the current row or column in the matrix
    maxIndex = 0  # the index of the max value
    for i in range((len(A))):  # calc the max value for each member on the main diagonal
        sum = 0  # the sum of the members in the current row in A
        for j in range(len(A)):  # go over the current row
            sum += abs(A[i][j])  # add the value of each member in the row
            if abs(A[i][j]) > max:  # search for the max value in the current row
                max = abs(A[i][j])
                maxIndex = j
        if (sum - max) <= max:  # if the max value in the row meets the condition of a dominant diagonal
            A = manualSwapCol(A, maxIndex,
                              i)  # swap between the columns of the current value on the main diagonal and the max value in that row
        else:  # look for the max value in the current column
            max = 0
            maxIndex = 0
            for j in range(len(A)):  # go over the current column
                # sum += abs(A[j][i])
                if abs(A[j][i]) > max:  # search for the max value in the current column
                    max = abs(A[j][i])
                    maxIndex = j
            if rowSum(A[j]) - max <= max:  # if the max value in the row meets the condition of a dominant diagonal
                A, b = manualSwapRow(A, b, i,
                                     maxIndex)  # swap between the rows of the current value on the main diagonal and the max value in that column
            else:
                print("ERROR - no dominant diagonal")  # A can't be changed into dominant diagonal matrix
                return None, None
    return A, b


def manualSwapRow(a, b, r1, r2):
    """
    manaul rows exchange (without e)
    :param a: The coefficient matrix
    :param b:  The column vector of constant terms
    :param r1: the first row to swap
    :param r2: the second row to swap
    :return: the matrix after the swap, The column vector of constant terms after swap
    """

    if r2 < len(a) and r1 < len(a):
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp
        if b is not None:  # if the result vector is not none swap him too
            temp = b[r1]
            b[r1] = b[r2]
            b[r2] = temp
    return a, b


def manualSwapCol(a, c1, c2):
    """
    :param a: The coefficient matrix
    :param c1: the first column to swap
    :param c2: the second column to swap
    :return: the matrix after the swap
    """
    if c2 < len(a) and c1 < len(a):
        for i in range(len(a)):
            temp = a[i][c1]
            a[i][c1] = a[i][c2]
            a[i][c2] = temp
    return a


def rowSum(line):
    """
    :param line: A list od numbers - line for the matrix
    :return: the sum of all the numbers in abs  in the list
    """
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


# end dominant part

# ############## main ###############
def driver():
    """
    main function for iterative isolation of variables method
    :return: prints results
    """
    A = [[1, 0, -1],
         [-0.5, 1, -0.25],
         [1, -0.5, 1]]

    b = [[0.2],
         [-1.425],
         [2]]

    epsilon = 10 ** (-4)

    flagDom = dominantDiagonal(A)  # check if there is a dominant diagonal

    # check
    if flagDom is False:
        copyA = copyMat(A)
        copyB = copyMat(b)
        copyA, copyB = createDominantDiagonal(copyA, copyB)  # change the matrix to be with dominant diagonal
        if (copyA is not None) and (copyB is not None):  # check if the return matrices are not none
            A = copyA
            b = copyB
            flagDom = True

        # end check
    LUdecomposition(A, b)
    iterativGuassSeidel(A, b, epsilon, flagDom)


driver()
