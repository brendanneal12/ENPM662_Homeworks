## Brendan Neal
## ENPM662 Homework 3 Question 1.1

##------------------------Importing Libraries-------------------------##
from sympy import *

init_printing(use_unicode=False, wrap_line=False)

##------------------------Variable Definitions------------------------##

#D-H Parameters
alpha, a, d, th = symbols('alpha a d theta')

#Symbolic Joint Variables for UR10
th1, th2, th3, th4, th5, th6 = symbols('th1 th2 th3 th4 th5 th6')


#Defining Parameters as found in Report in Arrays
alpha_array = [-pi/2, pi, pi, -pi/2, pi/2, 0]
theta_array = [th1, (pi/2)+th2, th3, (pi/2)+th4, th5, th6]
d_array = [128, 176, 163.9, 151.8, 115.7, 92.2] #mm
a_array = [ 0, -612.7, -571.6, 0, 0, 0] #mm

##------------Transformation Matrices Definitions--------------------##

Rz = Matrix([[cos(th), -sin(th), 0, 0],
             [sin(th),  cos(th), 0, 0],
             [      0,        0, 1, 0],
             [      0,        0, 0, 1]])

Tz = Matrix([[  1,  0,  0,  0],
             [  0,  1,  0,  0],
             [  0,  0,  1,  d],
             [  0,  0,  0,  1]])

Tx = Matrix([[  1,  0,  0,  a],
             [  0,  1,  0,  0],
             [  0,  0,  1,  0],
             [  0,  0,  0,  1]])

Rx = Matrix([[  1,        0,        0,  0],
             [  0,  cos(alpha), -sin(alpha),  0],
             [  0,  sin(alpha),  cos(alpha),  0],
             [  0,        0,        0,  1]])

##------------Computing Each Row's Transformation Matrix--------------##
T1 = Rz*Tz*Tx*Rx
T2 = Rz*Tz*Tx*Rx
T3 = Rz*Tz*Tx*Rx
T4 = Rz*Tz*Tx*Rx
T5 = Rz*Tz*Tx*Rx
T6= Rz*Tz*Tx*Rx

#Printing General Transformation Matris:
print('T_General:')
pprint(T6)
print('\n')

#Substituting Params into Each Transformation Matrix
T1 = T1.subs(alpha, alpha_array[0]).subs(th, theta_array[0]).subs(d, d_array[0]).subs(a, a_array[0])
T2 = T2.subs(alpha, alpha_array[1]).subs(th, theta_array[1]).subs(d, d_array[1]).subs(a, a_array[1])
T3 = T3.subs(alpha, alpha_array[2]).subs(th, theta_array[2]).subs(d, d_array[2]).subs(a, a_array[2])
T4 = T4.subs(alpha, alpha_array[3]).subs(th, theta_array[3]).subs(d, d_array[3]).subs(a, a_array[3])
T5 = T5.subs(alpha, alpha_array[4]).subs(th, theta_array[4]).subs(d, d_array[4]).subs(a, a_array[4])
T6 = T6.subs(alpha, alpha_array[5]).subs(th, theta_array[5]).subs(d, d_array[5]).subs(a, a_array[5])


##-----------------Printing Each Rows T Matrix------------------------##
print('DISPLAYING EACH INDIVIDUAL TRANSFORMATION MATRIX: \n')
print("T1 :")
pprint(T1)
print('\n')

print("T2 :")
pprint(T2)
print('\n')

print("T3 :")
pprint(T3)
print('\n')

print("T4 :")
pprint(T4)
print('\n')

print("T5 :")
pprint(T5)
print('\n')

print("T6 :")
pprint(T6)
print('\n')

##-----------------------Printing Final T Matrix------------------------##
T_Final = T1*T2*T3*T4*T5*T6
print('DISPLAYING FINAL T MATRIX: \n')
pprint(T_Final)
print('\n')

##-------------------------Geometric Validation--------------------------##
theta_v1_array = [pi/2, 0, 0, 0, 0, 0]
theta_v2_array = [0, pi/2, 0, 0, 0, 0]
theta_v3_array = [0, 0, pi/2, 0, 0, 0]
theta_v4_array = [0, 0, 0, pi/2, 0, 0]
theta_v5_array = [0, 0, 0, 0, pi/2, 0]

T_Val_1 = T_Final.subs(th1, theta_v1_array[0]).subs(th2, theta_v1_array[1]).subs(th3, theta_v1_array[2]).subs(th4, theta_v1_array[3]).subs(th5, theta_v1_array[4]).subs(th6, theta_v1_array[5])

T_Val_2 = T_Final.subs(th1, theta_v2_array[0]).subs(th2, theta_v2_array[1]).subs(th3, theta_v2_array[2]).subs(th4, theta_v2_array[3]).subs(th5, theta_v2_array[4]).subs(th6, theta_v2_array[5])

T_Val_3 = T_Final.subs(th1, theta_v3_array[0]).subs(th2, theta_v3_array[1]).subs(th3, theta_v3_array[2]).subs(th4, theta_v3_array[3]).subs(th5, theta_v3_array[4]).subs(th6, theta_v3_array[5])

T_Val_4 = T_Final.subs(th1, theta_v4_array[0]).subs(th2, theta_v4_array[1]).subs(th3, theta_v4_array[2]).subs(th4, theta_v4_array[3]).subs(th5, theta_v4_array[4]).subs(th6, theta_v4_array[5])

T_Val_5 = T_Final.subs(th1, theta_v5_array[0]).subs(th2, theta_v5_array[1]).subs(th3, theta_v5_array[2]).subs(th4, theta_v5_array[3]).subs(th5, theta_v5_array[4]).subs(th6, theta_v5_array[5])



##-----------------------Geometric Validation------------------------##
print('Theta 1 Rotated by 90:')
pprint(T_Val_1)
print('\n')

print('Theta 2 Rotated by 90:')
pprint(T_Val_2)
print('\n')

print('Theta 3 Rotated by 90:')
pprint(T_Val_3)
print('\n')

print('Theta 4 Rotated by 90:')
pprint(T_Val_4)
print('\n')

print('Theta 5 Rotated by 90:')
pprint(T_Val_5)
print('\n')









