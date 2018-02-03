#Prints the eigenvalue data that is output from ./QEQUICK total
#if gsl is not installed, use"plot_eigen_nogsl.py"
eigenvalue = open("eigenvalues.dat","r")
#eigenvector= open("eigenvector.dat",'r')
N = 100
# if for some reason your file does not have 100 eigenvalues, uncomment the lines below
#for i, l in enumerate(eigenvalue):
#    pass
#N = i + 1
EVALS = [[] for i in range(N)]
for i in range(0,N):
    Stringeigen = str(eigenvalue.readline())
    EVALS[i] = Stringeigen.split("\t")
for i in range(N):
    for j in range(len(EVALS[i])):
        EVALS[i][j] = float(EVALS[i][j])
eigenvalue.close()

import matplotlib.pyplot as plt

print("How many eigenvalues to plot")
nuser = input()
npr = int(nuser)
#for i in range(N):
for i in range(npr):
    if(i<3):
        nvals = [s+2 for s in range(len(EVALS[i])) ]
    else:
        nvals = [s+i for s in range(len(EVALS[i])) ]
    plt.plot(nvals,EVALS[i],'o', markersize=1)
plt.ylabel("Energy Eigenvalue")
plt.xlabel("Harmonic Oscillator Basis Number (dimension of the matrix)")
plottitle = "Plot of the first "+str(nuser) + " eigenvalues vs matrix size"
plt.title(plottitle)
plt.show()