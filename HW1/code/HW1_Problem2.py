## Brendan Neal
## ENPM662 Homework 1 Question 1.2

##------------------------Importing Libraries-------------------------##
import sympy
##--------------------Defining Symbolic Variables---------------------##

#Defining Time
t = sympy.symbols('t') #Time (s)

#Defining Constants
L2 = sympy.symbols('L2') #Length of Link 2
L3 = sympy.symbols('L3') #Length of Link 3
T2 = sympy.symbols('T2') #Theta 2

#Defining Paramaters as a Function of Time
L1 = sympy.Function('L1')(t) #Length of Link 1
T1 = sympy.Function('T1')(t) #Theta 1
T3 = sympy.Function('T3')(t) #Theta 2


##-------------Defining Expressions (Derived in Report)---------------##
X = L1*sympy.cos(T1) + L2*sympy.cos(T1 + T2) + L3*sympy.cos(T1 + T2 + T3) #X Positional Fwd Kinematics
Y = L1*sympy.sin(T1) + L2*sympy.sin(T1 + T2) + L3*sympy.sin(T1 + T2 + T3) #Y Positional Fwd Kinematics
Phi = T1 + T2 + T3

##-------------------------Taking Derivatives-------------------------##
X_Dot = sympy.diff(X,t) #X Velocity Fwd Kinematics
Y_Dot = sympy.diff(Y,t) #Y Velocity Fwd Kinematics
Phi_Dot = sympy.diff(Phi,t) #Phi Velocity Fwd Kinematics

##-------------------------Printing Results---------------------------##
print('Xdot: \n', sympy.pretty(X_Dot), '\n')
print('Ydot: \n', sympy.pretty(Y_Dot), '\n')
print('Phidot: \n', sympy.pretty(Phi_Dot), '\n')


##--------------------------------------------------------------------##
##---------------Deriving Inverse Velocity Kinematics-----------------##
##--------------------------------------------------------------------##

##--------------------Expanding Trig Expressions-----------------------##
X_Dot_Trig= sympy.expand_trig(X_Dot)
Y_Dot_Trig = sympy.expand_trig(Y_Dot)

##---------------------Simplifying Expressions-------------------------##
X_Dot_Expanded = X_Dot_Trig.expand()
Y_Dot_Expanded = Y_Dot_Trig.expand()
Phi_Dot_Expanded = Phi_Dot.expand()

#Collecting and Deriving Coefficients from the X_Dot Equatiion
T1_Dot_Collected_X = sympy.collect(X_Dot_Expanded, T1.diff())
T3_Dot_Collected_X = sympy.collect(X_Dot_Expanded, T3.diff())
L1_Dot_Collected_X = sympy.collect(X_Dot_Expanded, L1.diff())

Coeff_T1Dot_X = T1_Dot_Collected_X.coeff(T1.diff())
Coeff_T3Dot_X = T3_Dot_Collected_X.coeff(T3.diff())
Coeff_L1Dot_X = L1_Dot_Collected_X.coeff(L1.diff())



#Collecting and Deriving Coefficients from the Y_Dot Equatiion
T1_Dot_Collected_Y = sympy.collect(Y_Dot_Expanded, T1.diff())
T3_Dot_Collected_Y = sympy.collect(Y_Dot_Expanded, T3.diff())
L1_Dot_Collected_Y = sympy.collect(Y_Dot_Expanded, L1.diff())

Coeff_T1Dot_Y = T1_Dot_Collected_Y.coeff(T1.diff())
Coeff_T3Dot_Y = T3_Dot_Collected_Y.coeff(T3.diff())
Coeff_L1Dot_Y = L1_Dot_Collected_Y.coeff(L1.diff())

#Collecting and Deriving Coefficients from the Phi_Dot Equatiion
T1_Dot_Collected_Phi = sympy.collect(Phi_Dot_Expanded, T1.diff())
T3_Dot_Collected_Phi = sympy.collect(Phi_Dot_Expanded, T3.diff())
L1_Dot_Collected_Phi = sympy.collect(Phi_Dot_Expanded, L1.diff())


Coeff_T1Dot_Phi = T1_Dot_Collected_Phi.coeff(T1.diff())
Coeff_T3Dot_Phi = T3_Dot_Collected_Phi.coeff(T3.diff())
Coeff_L1Dot_Phi = L1_Dot_Collected_Phi.coeff(L1.diff())

##---------------Formulating Matrix for Inv Kinematics---------------##
A = sympy.Matrix([[Coeff_T1Dot_X,Coeff_T3Dot_X, Coeff_L1Dot_X],
                 [Coeff_T1Dot_Y, Coeff_T3Dot_Y, Coeff_L1Dot_Y],
                 [Coeff_T1Dot_Phi, Coeff_T3Dot_Phi, Coeff_L1Dot_Phi]])
A_inv = A.inv()

print("The inverse kinematic matrix is calculated but not displayed, please see the report to view the full parameters for the inverse kinematic matrix.")






















