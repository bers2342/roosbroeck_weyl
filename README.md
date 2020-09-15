# Numerical resolution of Van Roosbrock equations in a Weyl semimetal

The formating.py file is used for the plot part.

## Objective

This project is about computing the current, electrostatic potential and local chemical potential profile of different devices. One might consider the pesence of a voltage difference or a incident light beam. We consider the case of a 1D semiconductor and a 3D Weyl semimetal under strong magnetic field.
We consider only two Weyl nodes of opposite chirality.


## Code steps

### Quantities 

First we define useful quantities that we will use in the code. 
The main ones are given external or internal parameters that we input such as temperature *T*, magnetic field strength *B*, dieletric constant of the medium *\epsilon*, etc.
This is necessary to descritize in space our system. We focus on 1D system to descritize first because they are less expensive computationnaly. We choose an uniform spacing to start, as it is simpler.

### Van Roosbroeck equations

The next step is to write the Van Roosbroeck equations of the system. In the case without external light, if the Fermi level crosses only the zeroth Landau level, then only this level contributes.
We need three differential equations : one for the electrostatic potential (*F1*), one for the quasi-fermi level of the Left node (*F2_L*) and the last one for the Right node (*F2_R*).
Those notations are taken from the Chapter 50 of [this reference](https://www.taylorfrancis.com/books/e/9781315152318/chapters/10.4324/9781315152318-25)

### Solve the system of equations