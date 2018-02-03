#if you don't have gsl installed you can use this program to generate the eigen
#data from the file "hamiltonian.dat"
#this program assumes numpy (a common program) is installed

import numpy as np
import numpy.linalg as LA

eigenvalue = open("hamiltonian.dat","r")
#eigenvector= open("eigenvector.dat",'r')
N = 100
# if for some reason your file does not have 100 eigenvalues, uncomment the lines below
#for i, l in enumerate(eigenvalue):
#    pass
#N = i + 1
HAMILTONIAN = [[] for i in range(N)]
for i in range(0,N):
    Stringeigen = str(eigenvalue.readline())
    HAMILTONIAN[i] = Stringeigen.split("\t")
for i in range(N):
    for j in range(len(HAMILTONIAN[i])):
        HAMILTONIAN[i][j] = float(HAMILTONIAN[i][j])
eigenvalue.close()

Hnpy = np.array(HAMILTONIAN)
#print(Hnpy[0:3,0:3])
EVALS = [[] for i in range(N)]
for i in range(0,N):
    hsub = Hnpy[0:i+1 , 0:i+1]
    EVALS[i] = LA.eigvalsh(hsub)

#print(EVALS)
#print(EVALS[1])
#print(EVALS[N-1])

import matplotlib.pyplot as plt

print("How many eigenvalues to plot")
nuser = input()
npr = int(nuser)
#for i in range(N):
for i in range(npr):
    nvals = [ss+i for ss in range(N-i)]
    preval= [EVALS[ss+i][i] for ss in range(N-i)]
    #print(nvals)
    #print(preval)
    plt.plot(nvals,preval,'o', markersize=1)
plt.ylabel("Energy Eigenvalue")
plt.xlabel("Harmonic Oscillator Basis Number (dimension of the matrix)")
plottitle = "Plot of the first "+str(nuser) + " eigenvalues vs matrix size"
plt.title(plottitle)
plt.show()