#!/usr/bin/env python3

#Corresponding auhors : Didier Casane (didier.casane@egce.cnrs-gif.fr) and Maxime Policarpo (maxime.policarpo@egce.cnrs-gif.fr)


import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


T=79   #Number of genes
l=1045 #mean size of genes

dimension = (100000, (T+1))
bino_matrice = np.zeros(dimension)





for D in range(0,(T+1)):
	myvector=list()
	for time in range(1,1000001):
        	if(time % 10 == 0):
			X1=(((math.factorial(T))/((math.factorial(D)*(math.factorial(T-D))))))*((1-math.exp((-10**-8)*0.1278851*time*l))**D)*(math.exp((-10**-8)*0.1278851*time*l))**(T-D)
			myvector.append(X1)
	bino_matrice[:,D]=myvector





print(bino_matrice)

np.savetxt(("Matrice_loi_binomiale_lethops.csv"), bino_matrice, delimiter=",")

