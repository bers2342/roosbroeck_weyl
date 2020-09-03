# -*- coding: utf-8 -*-
import numpy as np
import sys
import matplotlib.pyplot as plt
import time
import formating as ft
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


###################
### FULL SOLVER ###
###################

'''

### ALL THESES STEPS ARE DESCRIBED IN THE FARREL PAPER ###


The first step is to declare all the variables, which are values of B, epsilon, etc.

Then we have to write the equations for psi, varphi^L and varphi^R. We use the central difference discretization, which converges quite well is this case.

Then we have to compute the Jacobian Matrix. In order to do that we need to compute the derivatives of the Van Roosboeck equations w.r.t. all the variables and insert them to the matrix

We have to write the loop in order to compute the acualized parameters in fonction of the old ones and the Jacobian matrix times the F-vector, which is a vector made of the Van Roosboeck equations

'''

# starting time
start = time.time()

# number of sites (Nsites = 501 is already a good converging result and the program has an exec time of 3 secs)
Nsites = 801

# uniform grid (in microns)
x = np.linspace(0., 10., Nsites)


# spacing between to points
h = x[1] - x[0]

# boundary subintervals
#w = 0.5*(x[1] + x[2]) - 0.5*(x[0] + x[1])
w = h


# magnetic field in [T]
B = 1.

# magnetic length in [microns]
ellB = 2.56*10**-2/np.sqrt(B) 

# Dielectric constant 
# We set epsilon >> 1 because we need a magnetic caracteristic length much bigger in order to have the code to converge
epsilon = 30.

# Fermi velocity in [mu m/s]
vF = 3.0 *10**11

#value of the tilt
tx = 0.
ty = 0.
tz = 0.
gamma = 1./np.sqrt(1. - tx**2 - ty**2) 

# relaxation time in [s]
tau_inter = 10**-12
tau_intra = 10**-12

# we define the constant : alpha_para = q^2 / (4*pi^2 * \epsilon_0 \epsilon \hbar v_F)
alpha_para = (1000/137)/(epsilon*np.pi)

# Landau gap in  [eV]
hvell = 7.65*10**-3*np.sqrt(B)

#print 'Landau gap = ', round(np.sqrt(2)*hvell*10**3, 3), 'meV'

print 'Magnetic caracteristic length = ', ellB/np.sqrt(alpha_para), 'microns'

print 'Mean free path = ', tau_inter*vF, 'microns'

print 'Spacing = ', h, 'microns'



### VAN ROOSBOECK EQUATIONS ###

###############################
### PARAMESWARAN LIKE PARAM ###
###############################

###########################################################
### THIS IS PARAMESWARAN PAPER PARAMETERS AND EQUATIONS ###
###########################################################

'''

not_qn = 1.
if not_qn == +1:
	print "Local neutrality ? : NO"
else :
	print "Local neutrality ? : YES"


# Discretized equation for psi(z) === also known as F1
def F1(psi_km1, psi_k, psi_kp1,phi0_L_km1, phi0_L_k, phi0_L_kp1, phi0_R_km1, phi0_R_k, phi0_R_kp1):
	return (psi_kp1 - 2.*psi_k + psi_km1)/h - alpha_para/ellB**2*w*(2.*psi_k - phi0_R_k - phi0_L_k)*not_qn

# Discretized equation for phi_L_0(z) === also known as F2_L
def F2_L(phi0_L_km1, phi0_L_k, phi0_L_kp1, psi_km1, psi_k, psi_kp1, phi0_R_km1, phi0_R_k, phi0_R_kp1):
	return tau_intra*vF*(phi0_L_kp1 - 2.*phi0_L_k + phi0_L_km1)/h**2 - 0.5*(phi0_L_kp1 - phi0_L_km1)/h - 1./(2.*vF*tau_inter)*(phi0_L_k - phi0_R_k)

# Discretized equation for phi_R_0(z) === also known as F2_R
def F2_R(phi0_R_km1, phi0_R_k, phi0_R_kp1, psi_km1, psi_k, psi_kp1, phi0_L_km1, phi0_L_k, phi0_L_kp1):
	return tau_intra*vF*(phi0_R_kp1 - 2.*phi0_R_k + phi0_R_km1)/h**2 + 0.5*(phi0_R_kp1 - phi0_R_km1)/h + 1./(2.*vF*tau_inter)*(phi0_L_k - phi0_R_k)



### COMPUTATIONS OF DERIVATIVES FOR THE JACOBIAN MATRIX ###

### F1 ###


def dF1_dpsi_k(psi_k):
	return -2./h - 2.*w*alpha_para/ellB**2*not_qn

def dF1_dpsi_kpm1(psi_k):
	return 1./h

def dF1_dphi0L_k(psi_k):
	return w*alpha_para/ellB**2*not_qn

def dF1_dphi0R_k(psi_k):
	return w*alpha_para/ellB**2*not_qn 


def dF1_dphi0L_kp1(psi_k):
	return 0.

def dF1_dphi0L_km1(psi_k):
	return 0.


def dF1_dphi0R_kp1(psi_k):
	return 0.

def dF1_dphi0R_km1(psi_k):
	return 0.

### F2L ### (---> first without E_ext)


def dF2L_dphi0L_k(phi0_L_k):
	return -2.*tau_intra*vF/h**2 - 1./(2.*vF*tau_inter)

def dF2L_dphi0L_kp1(phi0_L_kp1):
	return tau_intra*vF/h**2 - 0.5/h

def dF2L_dphi0L_km1(phi0_L_km1):
	return tau_intra*vF/h**2 + 0.5/h



def dF2L_dphi0R_k(phi0_R_k):
	return 1./(2.*vF*tau_inter)


def dF2L_dphi0R_kp1(phi0_R_k):
	return 0.

def dF2L_dphi0R_km1(phi0_R_k):
	return 0.

### F2R ###


def dF2R_dphi0R_k(phi0_L_k):
	return -2.*tau_intra*vF/h**2 - 1./(2.*vF*tau_inter)

def dF2R_dphi0R_kp1(phi0_L_kp1):
	return tau_intra*vF/h**2 + 0.5/h

def dF2R_dphi0R_km1(phi0_L_km1):
	return tau_intra*vF/h**2 - 0.5/h 



def dF2R_dphi0L_k(phi0_R_k):
	return 1./(2.*vF*tau_inter)


def dF2R_dphi0L_kp1(phi0_L_kp1):
	return 0.

def dF2R_dphi0L_km1(phi0_L_kp1):
	return 0.


'''

###############################
### VAN ROOSBOECK EQUATIONS ###
###############################

###############################
### PARAMESWARAN LIKE PARAM ###
###############################


##############################################
########## QUASINEUTRAL CASE OR NOT ##########
### not_qn = +1 for NOT QN and  = 0 for QN ###
##############################################

not_qn = +1
if not_qn == +1:
	print "Local neutrality ? : NO"
else :
	print "Local neutrality ? : YES"


# Discretized equation for psi(z) === also known as F1
def F1(psi_km1, psi_k, psi_kp1, phi0_L_km1, phi0_L_k, phi0_L_kp1, phi0_R_km1, phi0_R_k, phi0_R_kp1):
	return (psi_kp1 - 2.*psi_k + psi_km1)/h - w*alpha_para/ellB**2*((2.*psi_k - phi0_R_k - phi0_L_k)/(1./gamma - tz) - 0.5*vF*tau_intra*(phi0_L_kp1 - phi0_L_km1)/h + 0.5*vF*tau_intra*(phi0_R_kp1 - phi0_R_km1)/h)*not_qn

# Discretized equation for phi_L_0(z) === also known as F2_L
def F2_L(phi0_L_km1, phi0_L_k, phi0_L_kp1, psi_km1, psi_k, psi_kp1, phi0_R_km1, phi0_R_k, phi0_R_kp1):
	return tau_intra*vF*(phi0_L_kp1 - 2.*phi0_L_k + phi0_L_km1)/h**2 + 0.5*(phi0_L_kp1 - phi0_L_km1)/h - 1./(2.*vF*tau_inter)*((phi0_L_k - phi0_R_k)/(1./gamma - tz) + vF*tau_intra*0.5*(phi0_L_kp1 - phi0_L_km1)/h + vF*tau_intra*0.5*(phi0_R_kp1 - phi0_R_km1)/h)

# Discretized equation for phi_R_0(z) === also known as F2_R
def F2_R(phi0_R_km1, phi0_R_k, phi0_R_kp1, psi_km1, psi_k, psi_kp1, phi0_L_km1, phi0_L_k, phi0_L_kp1):
	return tau_intra*vF*(phi0_R_kp1 - 2.*phi0_R_k + phi0_R_km1)/h**2 - 0.5*(phi0_R_kp1 - phi0_R_km1)/h + 1./(2.*vF*tau_inter)*((phi0_L_k - phi0_R_k)/(1./gamma - tz) + vF*tau_intra*0.5*(phi0_L_kp1 - phi0_L_km1)/h + vF*tau_intra*0.5*(phi0_R_kp1 - phi0_R_km1)/h)



### COMPUTATIONS OF DERIVATIVES FOR THE JACOBIAN MATRIX ###
################### DERIVATIVES ###########################


### F1 ###

def dF1_dpsi_k(psi_k):
	return -2./h - 2./(1./gamma - tz)*w*alpha_para/ellB**2*not_qn

def dF1_dpsi_kpm1(psi_k):
	return 1./h

def dF1_dphi0L_k(psi_k):
	return 1./(1./gamma - tz)*w*alpha_para/ellB**2*not_qn

def dF1_dphi0R_k(psi_k):
	return 1./(1./gamma - tz)*w*alpha_para/ellB**2*not_qn



def dF1_dphi0L_kp1(psi_k):
	return - 0.5*alpha_para/ellB**2*vF*tau_intra*not_qn

def dF1_dphi0L_km1(psi_k):
	return + 0.5*alpha_para/ellB**2*vF*tau_intra*not_qn


def dF1_dphi0R_kp1(psi_k):
	return + 0.5*alpha_para/ellB**2*vF*tau_intra*not_qn

def dF1_dphi0R_km1(psi_k):
	return - 0.5*alpha_para/ellB**2*vF*tau_intra*not_qn


### F2L ### 


def dF2L_dphi0L_k(phi0_L_k):
	return -2.*tau_intra*vF/h**2 - 1./(2.*vF*tau_inter*(1./gamma - tz))

def dF2L_dphi0L_kp1(phi0_L_kp1):
	return tau_intra*vF/h**2 + 0.5/h - 0.5/(2.*vF*tau_inter)*vF*tau_intra/h

def dF2L_dphi0L_km1(phi0_L_km1):
	return tau_intra*vF/h**2 - 0.5/h + 0.5/(2.*vF*tau_inter)*vF*tau_intra/h




def dF2L_dphi0R_k(phi0_R_k):
	return 1./(2.*vF*tau_inter*(1./gamma - tz))

def dF2L_dphi0R_kp1(phi0_R_k):
	return - 0.5/(2.*vF*tau_inter)*(vF*tau_intra/h)

def dF2L_dphi0R_km1(phi0_R_k):
	return + 0.5/(2.*vF*tau_inter)*(vF*tau_intra/h)


### F2R ###


def dF2R_dphi0R_k(phi0_L_k):
	return -2.*tau_intra*vF/h**2 - 1./(2.*vF*tau_inter*(1./gamma - tz))

def dF2R_dphi0R_kp1(phi0_L_kp1):
	return tau_intra*vF/h**2 - 0.5/h + 0.5/(2.*vF*tau_inter)*vF*tau_intra/h

def dF2R_dphi0R_km1(phi0_L_km1):
	return tau_intra*vF/h**2 + 0.5/h - 0.5/(2.*vF*tau_inter)*vF*tau_intra/h




def dF2R_dphi0L_k(phi0_L_k):
	return 1./(2.*vF*tau_inter*(1./gamma - tz))

def dF2R_dphi0L_kp1(phi0_L_kp1):
	return + 0.5/(2.*vF*tau_inter)*(vF*tau_intra/h)

def dF2R_dphi0L_km1(phi0_L_kp1):
	return - 0.5/(2.*vF*tau_inter)*(vF*tau_intra/h)



### FULL SYSTEM JACOBIAN MATRIX ###

# boundary conditions for VARPHI in [V] or [eV]

# step of the voltage difference
dU = 0.1

# Initilization of the edges
psi_1 = 0.
psi_N = 0.

phi0_L_1 = 0.
phi0_L_N = 0.

phi0_R_1 = 0.
phi0_R_N = 0.


# initial values in [V] or [eV]
psi_0 = np.zeros(Nsites)
phi0_L_0 = np.zeros(Nsites)
phi0_R_0 = np.zeros(Nsites)

#psi_0 = (psi_N-psi_1)*x/x[-1] + psi_1
#phi0_L_0 = (psi_N-psi_1)*x/x[-1] + psi_1
#phi0_R_0 = (psi_N-psi_1)*x/x[-1] + psi_1


### JACOBIAN MATRIX ###
def full_Jacobian_WSM(psi, phi0_L, phi0_R):
	jac_mat = np.zeros((3*Nsites, 3*Nsites))
	for i in range(Nsites):
		for j in range(Nsites):
			if i==j:
				jac_mat[i][j] = dF1_dpsi_k(psi[i])
				jac_mat[i][j + Nsites] = dF1_dphi0L_k(psi[i])
				jac_mat[i][j + 2*Nsites] = dF1_dphi0R_k(psi[i])
				
				jac_mat[i + Nsites][j + Nsites] = dF2L_dphi0L_k(phi0_L[i])
				jac_mat[i + Nsites][j + 2*Nsites] = dF2L_dphi0R_k(phi0_L[i])

				jac_mat[i + 2*Nsites][j + 2*Nsites] = dF2R_dphi0R_k(phi0_R[i])
				jac_mat[i + 2*Nsites][j + Nsites] = dF2R_dphi0L_k(phi0_L[i])

			if abs(i-j) == 1:
				jac_mat[i][j] = dF1_dpsi_kpm1(psi[i])

			if j == i+1:
				jac_mat[i][j + Nsites] = dF1_dphi0L_kp1(psi[i])
				jac_mat[i][j + 2*Nsites] = dF1_dphi0R_kp1(psi[i])

				jac_mat[i + Nsites][j + 2*Nsites] = dF2L_dphi0R_kp1(phi0_R[i])
				jac_mat[i + Nsites][j + Nsites] = dF2L_dphi0L_kp1(phi0_L[i])

				jac_mat[i + 2*Nsites][j + Nsites] = dF2R_dphi0L_kp1(phi0_L[i])
				jac_mat[i + 2*Nsites][j + 2*Nsites] = dF2R_dphi0R_kp1(phi0_R[i])

			if j == i-1:
				jac_mat[i][j + Nsites] = dF1_dphi0L_km1(psi[i])
				jac_mat[i][j + 2*Nsites] = dF1_dphi0R_km1(psi[i])

				jac_mat[i + Nsites][j + 2*Nsites] = dF2L_dphi0R_km1(phi0_R[i])
				jac_mat[i + Nsites][j + Nsites] = dF2L_dphi0L_km1(phi0_L[i])

				jac_mat[i + 2*Nsites][j + Nsites] = dF2R_dphi0L_km1(phi0_L[i])
				jac_mat[i + 2*Nsites][j + 2*Nsites] = dF2R_dphi0R_km1(phi0_R[i])		
	return jac_mat

#print 'Jacobian Matrix = ', full_Jacobian_WSM(psi0, phi0_L_0, phi0_R_0, phi1_L_0, phi1_R_0)


############## FULL SYSTEM F VECTOR #############
### WE HAVE TO SPCIFY THE BOUNDARY CONDITIONS ###

def full_Fvector_WSM(psi, phi0_L, phi0_R):
	f_vec = np.zeros(3*Nsites)
	for i in range(Nsites):
		if i==0:
			f_vec[i] = F1(psi_1, psi[i], psi[i+1], phi0_L_1, phi0_L[i], phi0_L[i+1], phi0_R_1, phi0_R[i], phi0_R[i+1])
			f_vec[i + Nsites] = F2_L(phi0_L_1, phi0_L[i], phi0_L[i+1], psi_1, psi[i], psi[i+1], phi0_R_1, phi0_R[i], phi0_R[i+1])
			f_vec[i + 2*Nsites] = F2_R(phi0_R_1, phi0_R[i], phi0_R[i+1], psi_1, psi[i], psi[i+1], phi0_L_1, phi0_L[i], phi0_L[i+1])
		elif i==Nsites-1:
			f_vec[i] = F1(psi[i-1], psi[i], psi_N, phi0_L[i-1], phi0_L[i], phi0_L_N, phi0_R[i-1], phi0_R[i], phi0_R_N)
			f_vec[i + Nsites] = F2_L(phi0_L[i-1], phi0_L[i], phi0_L_N, psi[i-1], psi[i], psi_N, phi0_R[i-1], phi0_R[i], phi0_R_N)
			f_vec[i + 2*Nsites] = F2_R(phi0_R[i-1], phi0_R[i], phi0_R_N, psi[i-1], psi[i], psi_N, phi0_L[i-1], phi0_L[i], phi0_L_N)
		else:
			f_vec[i] = F1(psi[i-1], psi[i], psi[i+1], phi0_L[i-1], phi0_L[i], phi0_L[i+1], phi0_R[i-1], phi0_R[i], phi0_R[i+1])
			f_vec[i + Nsites] = F2_L(phi0_L[i-1], phi0_L[i], phi0_L[i+1], psi[i-1], psi[i], psi[i+1], phi0_R[i-1], phi0_R[i], phi0_R[i+1])
			f_vec[i + 2*Nsites] = F2_R(phi0_R[i-1], phi0_R[i], phi0_R[i+1], psi[i-1], psi[i], psi[i+1], phi0_L[i-1], phi0_L[i], phi0_L[i+1])
	return f_vec


### NEWTON'S LOOP ### 


# Loop iterations
# We define the initialization vector from the initial vector values
vec = []
# The vector is the concatenation of psi, varphi^L and varphi^R
vec0 = np.concatenate((psi_0, phi0_L_0, phi0_R_0), axis=None)
vec.append(vec0)


# These temp variables will be the initialized values in order to compute the first iteration of the matrix
vec_temp = vec0
psi_temp, phi0_L_temp, phi0_R_temp = psi_0, phi0_L_0, phi0_R_0

### CONVERGENCE LOOP ###
for i in range(10):
	# This is the update step where psi_update = psi_old - Jacobian * F (same for varphi)
	vec_iter = vec_temp - np.dot(np.linalg.inv(full_Jacobian_WSM(psi_temp, phi0_L_temp, phi0_R_temp)),full_Fvector_WSM(psi_temp, phi0_L_temp, phi0_R_temp))
	vec.append(vec_iter)
	# WE HAVE TO UPDATE THE BOUNDARY VOLTAGES
	if psi_1 < 5.*dU :
		psi_1    += dU
		phi0_L_1 += dU
		phi0_R_1 += dU
		psi_N    -= dU
		phi0_L_N -= dU
		phi0_R_N -= dU
	# WE SAVE THE NEW VALUES
	psi_temp, phi0_L_temp, phi0_R_temp = vec_iter[:Nsites], vec_iter[Nsites:2*Nsites], vec_iter[2*Nsites:]
	# THE UPDATED VECTOR BECOMES THE NEW INITIALZED VECTOR AND THE LOOP GOES ON.
	vec_temp = vec_iter

### NOTE ###
'''
It is very important to update the boundary values in the loop in order for the Jacobian and the F vector to have the GOOD boundary values for the BIAS
'''

# current in units of q^2/(4 pi^2 ellB^2 hbar)
def j0chi(psi, phi0, chi):
	j0chi = np.zeros(Nsites)
	for i in range(Nsites):
		if chi == 1.:
			if i == Nsites-1:
				j0chi[i] = (psi[i] - phi0[i]) - 0.5*vF*tau_intra*(phi0_L_N - phi0[i-1])/h
			elif i==0:
				j0chi[i] = (psi[i] - phi0[i]) - 0.5*vF*tau_intra*(phi0[i+1] - phi0_L_1)/h
			else:
				j0chi[i] = (psi[i] - phi0[i]) - 0.5*vF*tau_intra*(phi0[i+1] - phi0[i-1])/h
		if chi == -1. :
			if i== Nsites-1:
				j0chi[i] = -(psi[i] - phi0[i]) - 0.5*vF*tau_intra*(phi0_R_N - phi0[i-1])/h
			elif i==0:
				j0chi[i] = -(psi[i] - phi0[i]) - 0.5*vF*tau_intra*(phi0[i+1] - phi0_R_1)/h
			else:
				j0chi[i] = -(psi[i] - phi0[i]) - 0.5*vF*tau_intra*(phi0[i+1] - phi0[i-1])/h
	return j0chi



print 'time = ', time.time() - start, 's'

'''
#plt.plot(x, vec[0][:Nsites], ':',   linewidth=3.0, color='black', label=r"$ start $")
plt.plot(x, vec[1][:Nsites], '-.',  linewidth=3.5, color='gold', label=r" Iter 0 ")
plt.plot(x, vec[2][:Nsites], ':',   linewidth=3.5, color='green', label=r" Iter 1 ")
plt.plot(x, vec[3][:Nsites], '-',   linewidth=3.5, color='blue', label=r" Iter 2 ")
plt.plot(x, vec[-1][:Nsites], '--', linewidth=3.5, color='red', label=r" Iter "+str(len(vec)-1)+" ")
plt.plot(x, (psi_N-psi_1)*x/x[Nsites-1] + psi_1, '-', linewidth=2.0, color='black', label=r" Res. ana. ")
plt.xlim(x[0], x[Nsites-1])
plt.yticks(np.arange(-0.5, 0.6, 0.1))
plt.xlabel(r"$ z [\mu m] $")
plt.ylabel(r"$ \psi $[mV]")
plt.legend()
plt.show()

#plt.plot(x, vec[0][Nsites:2*Nsites], ':',   linewidth=3.0, color='black', label=r"$ 0^{th} $ iter ")
plt.plot(x[::15], vec[1][Nsites:2*Nsites][::15], 's',  ms=5.0, color='gold', label=r" Iter 0")
plt.plot(x[::15], vec[2][Nsites:2*Nsites][::15], 'o',   ms=5.0, color='green', label=r" Iter 1")
plt.plot(x[::15], vec[3][Nsites:2*Nsites][::15], 'v',   ms=5.0, color='blue', label=r" Iter 2")
plt.plot(x[::15], vec[-1][Nsites:2*Nsites][::15], '^', ms=5.0, color='red', label=r" Iter "+str(len(vec)-1)+" ")


#plt.xlim(x[0], x[Nsites-1])
#plt.xlabel(r"$ z [\mu m] $")
#plt.ylabel(r"$ \varphi_0^L [V]$")
#plt.legend()
#plt.show()


#plt.plot(x, vec[0][2*Nsites:], ':',   linewidth=3.0, color='black', label=r"$ 0^{th} $ iter ")
#plt.plot(x, vec[1][2*Nsites:], '-.',  linewidth=3.0, color='gold')
#plt.plot(x, vec[2][2*Nsites:], ':',   linewidth=3.0, color='green')
#plt.plot(x, vec[3][2*Nsites:], '-',   linewidth=3.0, color='blue')
#plt.plot(x, vec[-1][2*Nsites:], '--', linewidth=3.0, color='red')
plt.plot(x, (psi_N-psi_1)*x/x[Nsites-1] + psi_1, '-', linewidth=2.0, color='black', label=r" Res. ana.")

plt.xlim(x[0], x[Nsites-1])
plt.yticks(np.arange(-0.5, 0.6, 0.1))
plt.xlabel(r"$ z [\mu m] $")
plt.ylabel(r"$ \varphi_0 $[mV]")
plt.legend()
plt.show()




plt.plot(x, j0chi(vec[-1][:Nsites], vec[-1][Nsites:2*Nsites], 1.), '-.', linewidth=4.0, color='red', label=r"$ j_0^{+} (z) $")
plt.plot(x, j0chi(vec[-1][:Nsites], vec[-1][2*Nsites:], -1.), ':', linewidth=4.0, color='blue', label=r"$ j_0^{-} (z) $")

plt.plot(x, j0chi(vec[-1][:Nsites], vec[-1][Nsites:2*Nsites], 1.) + j0chi(vec[-1][:Nsites], vec[-1][2*Nsites:], -1.), '-', linewidth=4.0, color='black', label=r"$ j_0^{c} (z) $")
plt.plot(x, j0chi(vec[-1][:Nsites], vec[-1][Nsites:2*Nsites], 1.) - j0chi(vec[-1][:Nsites], vec[-1][2*Nsites:], -1.), '--', linewidth=4.0, color='gold', label=r"$ j_0^{ 5 } (z) $")


plt.xlim(x[0], x[Nsites-1])
plt.xlabel(r"$ z [\mu m] $")
plt.ylabel(r"$ j_0 [q^2/(4 \pi^2 \ell_B^2 \hbar)]$")
plt.legend(loc='upper right')
plt.show()
'''


fontsize=25
pdftype=3

#figures_path = "psi_verfied_v2.pdf"
#figures_path = "varphi_verfied_v2.pdf"
figures_path = "current_verfied_v2.pdf"

xmin,xmax=x[0], x[-1]
ymin,ymax= -0.003, 0.063

space_x=2.
space_y=0.01

xlabel=r"$z [ \mu m]$"
ylabel=r"$j_0 [q^2 / (4 \pi^2 \ell_B^2 \hbar)]$"
title=""

## Create fig and axes ################################

ft.step1(fontsize,pdftype)
fig, axes = plt.subplots(1,1,figsize=(8.5, 6)) #create canvas
fig.subplots_adjust(left=0.17, right=0.95, bottom=0.15, top=0.95) #adjust fig on canvas

## Mise en forme des axes #################

ft.set_ax(axes,xmin,xmax,ymin,ymax,xlabel,ylabel,title)
ft.set_ticks(axes,space_x,space_y)

## Let's plot ###############################

#axes.plot(x, vec[1][:Nsites], '-.',  linewidth=4., color='gold', label=r" Iter 0 ")
#axes.plot(x, vec[2][:Nsites], ':',   linewidth=4., color='green', label=r" Iter 1 ")
#axes.plot(x, vec[3][:Nsites], '-',   linewidth=4., color='blue', label=r" Iter 2 ")
#axes.plot(x, vec[-1][:Nsites], '--', linewidth=4., color='red', label=r" Iter "+str(len(vec)-1)+" ")
#axes.plot(x, (psi_N-psi_1)*x/x[Nsites-1] + psi_1, '-', linewidth=2.0, color='black', label=r" Res. ana. ")


#axes.plot(x[::15], vec[1][Nsites:2*Nsites][::15], 's',  ms=7.0, color='gold', label=r" Iter 0")
#axes.plot(x[::15], vec[2][Nsites:2*Nsites][::15], 'o',   ms=7.0, color='green', label=r" Iter 1")
#axes.plot(x[::15], vec[3][Nsites:2*Nsites][::15], 'v',   ms=7.0, color='blue', label=r" Iter 2")
#axes.plot(x[::15], vec[-1][Nsites:2*Nsites][::15], '^', ms=7.0, color='red', label=r" Iter "+str(len(vec)-1)+" ")
#axes.plot(x, (psi_N-psi_1)*x/x[Nsites-1] + psi_1, '-', linewidth=2.0, color='black', label=r" Res. ana.")

axes.plot(x, j0chi(vec[-1][:Nsites], vec[-1][Nsites:2*Nsites], 1.), '-.', linewidth=4.5, color='red', label=r"$ j_0^{+} (z) $")
axes.plot(x, j0chi(vec[-1][:Nsites], vec[-1][2*Nsites:], -1.), ':', linewidth=4.5, color='blue', label=r"$ j_0^{-} (z) $")

axes.plot(x, j0chi(vec[-1][:Nsites], vec[-1][Nsites:2*Nsites], 1.) + j0chi(vec[-1][:Nsites], vec[-1][2*Nsites:], -1.), '-', linewidth=4.5, color='black', label=r"$ j_0^{c} (z) $")
axes.plot(x, j0chi(vec[-1][:Nsites], vec[-1][Nsites:2*Nsites], 1.) - j0chi(vec[-1][:Nsites], vec[-1][2*Nsites:], -1.), '--', linewidth=4.5, color='gold', label=r"$ j_0^{5 } (z) $")


plt.legend(loc=1,fontsize=18)

s=20


plt.show()
fig.savefig(figures_path)