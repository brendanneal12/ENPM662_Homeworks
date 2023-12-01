# Brendan Neal
# ENPM662 Homework 5

'''Please note I am reusing some code from HW3 and HW4'''

## ------------------------Importing Libraries-------------------------##
from sympy import *
from matplotlib import pyplot as plt
import numpy as np

init_printing(use_unicode=False, wrap_line=False)

## ------------------------Variable Definitions------------------------##

# D-H Parameters
alpha, a, d, th = symbols('alpha a d theta')

# Symbolic Joint Variables for UR10
th1, th2, th3, th4, th5, th6 = symbols('th1 th2 th3 th4 th5 th6')


# Defining Parameters as found in Report in Arrays
alpha_array = [-pi/2, pi, pi, -pi/2, pi/2, 0]
theta_array = [th1, (pi/2)+th2, th3, (pi/2)+th4, th5, th6]
d_array = [.128, .176, .1639, .1518, .1157, .1922]  # m
a_array = [0, -.6127, -.5716, 0, 0, 0]  # m

# Mass Information (kg)
m1 = 7.1
m2 = 12.7
m3 = 4.27
m4 = 2
m5 = 2
m6 = 0.365

# Gravity Info
g = 9.81  # m/s^2

## ------------Transformation Matrices Definitions--------------------##

Rz = Matrix([[cos(th), -sin(th), 0, 0],
             [sin(th),  cos(th), 0, 0],
             [0,        0, 1, 0],
             [0,        0, 0, 1]])

Tz = Matrix([[1,  0,  0,  0],
             [0,  1,  0,  0],
             [0,  0,  1,  d],
             [0,  0,  0,  1]])

Tx = Matrix([[1,  0,  0,  a],
             [0,  1,  0,  0],
             [0,  0,  1,  0],
             [0,  0,  0,  1]])

Rx = Matrix([[1,        0,        0,  0],
             [0,  cos(alpha), -sin(alpha),  0],
             [0,  sin(alpha),  cos(alpha),  0],
             [0,        0,        0,  1]])

## ------------Computing Each Row's Transformation Matrix--------------##
T1_gen = Rz*Tz*Tx*Rx
T2_gen = Rz*Tz*Tx*Rx
T3_gen = Rz*Tz*Tx*Rx
T4_gen = Rz*Tz*Tx*Rx
T5_gen = Rz*Tz*Tx*Rx
T6_gen = Rz*Tz*Tx*Rx

# Substituting Params into Each Transformation Matrix
T1 = T1_gen.subs(alpha, alpha_array[0]).subs(
    th, theta_array[0]).subs(d, d_array[0]).subs(a, a_array[0])
T2 = T2_gen.subs(alpha, alpha_array[1]).subs(
    th, theta_array[1]).subs(d, d_array[1]).subs(a, a_array[1])
T3 = T3_gen.subs(alpha, alpha_array[2]).subs(
    th, theta_array[2]).subs(d, d_array[2]).subs(a, a_array[2])
T4 = T4_gen.subs(alpha, alpha_array[3]).subs(
    th, theta_array[3]).subs(d, d_array[3]).subs(a, a_array[3])
T5 = T5_gen.subs(alpha, alpha_array[4]).subs(
    th, theta_array[4]).subs(d, d_array[4]).subs(a, a_array[4])
T6 = T6_gen.subs(alpha, alpha_array[5]).subs(
    th, theta_array[5]).subs(d, d_array[5]).subs(a, a_array[5])

## ----------Defining Transformation Matrices wrt Zero Frame-----------##
H0_1 = T1
H0_2 = T1*T2
H0_3 = T1*T2*T3
H0_4 = T1*T2*T3*T4
H0_5 = T1*T2*T3*T4*T5
H0_6 = T1*T2*T3*T4*T5*T6

## ----------------------Extracting End Effector Position from H0_6------------------------##
P = Matrix([[H0_6[3]], [H0_6[7]], [H0_6[11]]])

## -----------Computing General Form of Partial Derivatives------------##

Par1 = diff(P, th1)
Par2 = diff(P, th2)
Par3 = diff(P, th3)
Par4 = diff(P, th4)
Par5 = diff(P, th5)
Par6 = diff(P, th6)

## -----------------Extracting Z from each H Matrix--------------------##
Z0_1 = Matrix([[H0_1[2]], [H0_1[6]], [H0_1[10]]])
Z0_2 = Matrix([[H0_2[2]], [H0_2[6]], [H0_2[10]]])
Z0_3 = Matrix([[H0_3[2]], [H0_3[6]], [H0_3[10]]])
Z0_4 = Matrix([[H0_4[2]], [H0_4[6]], [H0_4[10]]])
Z0_5 = Matrix([[H0_5[2]], [H0_5[6]], [H0_5[10]]])
Z0_6 = Matrix([[H0_6[2]], [H0_6[6]], [H0_6[10]]])

## -------------------Forming Jacobian Matrix---------------------##

# Components
J1 = Matrix([[Par1], [Z0_1]])
J2 = Matrix([[Par2], [Z0_2]])
J3 = Matrix([[Par3], [Z0_3]])
J4 = Matrix([[Par4], [Z0_4]])
J5 = Matrix([[Par5], [Z0_5]])
J6 = Matrix([[Par6], [Z0_6]])

# Full Jacobian
J = Matrix([[J1, J2, J3, J4, J5, J6]])

## ------------------Calculating Gravity Matrix-------------------##

# Calculating Center of Mass Information
d_array_new_1 = [.064, .176, .1639, .1518, .1157, .1922]  # m
a_array_new_1 = [0, -.6127, -.5716, 0, 0, 0]  # m

d_array_new_2 = [.128, .088, .1639, .1518, .1157, .1922]  # m
a_array_new_2 = [0, -.30635, -.5716, 0, 0, 0]  # m

d_array_new_3 = [.128, .176, .08195, .1518, .1157, .1922]  # m
a_array_new_3 = [0, -.6127, -.2858, 0, 0, 0]  # m

d_array_new_4 = [.128, .176, .1639, .0759, .1157, .1922]  # m
a_array_new_4 = [0, -.6127, -.5716, 0, 0, 0]  # m

d_array_new_5 = [.128, .176, .1639, .1518, .05785, .1922]  # m
a_array_new_5 = [0, -.6127, -.5716, 0, 0, 0]  # m

d_array_new_6 = [.128, .176, .1639, .1518, .1157, .0961]  # m
a_array_new_6 = [0, -.6127, -.5716, 0, 0, 0]  # m

T1_new = T1_gen.subs(alpha, alpha_array[0]).subs(
    th, theta_array[0]).subs(d, d_array_new_1[0]).subs(a, a_array_new_1[0])
T2_new = T2_gen.subs(alpha, alpha_array[1]).subs(
    th, theta_array[1]).subs(d, d_array_new_2[1]).subs(a, a_array_new_2[1])
T3_new = T3_gen.subs(alpha, alpha_array[2]).subs(
    th, theta_array[2]).subs(d, d_array_new_3[2]).subs(a, a_array_new_3[2])
T4_new = T4_gen.subs(alpha, alpha_array[3]).subs(
    th, theta_array[3]).subs(d, d_array_new_4[3]).subs(a, a_array_new_4[3])
T5_new = T5_gen.subs(alpha, alpha_array[4]).subs(
    th, theta_array[4]).subs(d, d_array_new_5[4]).subs(a, a_array_new_5[4])
T6_new = T6_gen.subs(alpha, alpha_array[5]).subs(
    th, theta_array[5]).subs(d, d_array_new_6[5]).subs(a, a_array_new_6[5])

H0_1_new = T1_new
H0_2_new = T1*T2_new
H0_3_new = T1*T2*T3_new
H0_4_new = T1*T2*T3*T4_new
H0_5_new = T1*T2*T3*T4*T5_new
H0_6_new = T1*T2*T3*T4*T5*T6_new

# Potential Energy Components
P1 = m1 * g * H0_1_new[11]
P2 = m2 * g * H0_2_new[11]
P3 = m3 * g * H0_3_new[11]
P4 = m4 * g * H0_4_new[11]
P5 = m5 * g * H0_5_new[11]
P6 = m6 * g * H0_6_new[11]

# Total Potential Energy

Ptot = P1 + P2 + P3 + P4 + P5 + P6

# Gravity Matrix Components
G1 = diff(Ptot, th1)
G2 = diff(Ptot, th2)
G3 = diff(Ptot, th3)
G4 = diff(Ptot, th4)
G5 = diff(Ptot, th5)
G6 = diff(Ptot, th6)

# Full Gravity Matrix
Gravity = Matrix([[G1], [G2], [G3], [G4], [G5], [G6]])
print("Gravity Matrix:")
pprint(Gravity)

## -------------------Plotting the Circle Code---------------------##

r = 0.1  # m

# Defining my initial joint angles.
q1 = 0.0002
q2 = 0.0001
q3 = -0.0001
q4 = 0.0001
q5 = 0.0004
q6 = 0.00001

# Force Matrix
F = Matrix([[0], [-5], [0], [0], [0], [0]])

# Torque Arrays for Plotting
Torque1 = []
Torque2 = []
Torque3 = []
Torque4 = []
Torque5 = []
Torque6 = []


time = np.linspace(0, 200, num=400)
dt = 200 / 400
for t in time:
    # Define components of X_Dot
    Vx = ((2*pi*r)/200)*cos((2*pi*t)/200)
    Vy = 0.0
    Vz = -((2*pi*r)/200)*sin((2*pi*t)/200)
    omega_x = 0.0
    omega_y = 0.0
    omega_z = 0.0

    # Formulate X_Dot
    X_dot = Matrix([[Vx], [Vy], [Vz], [omega_x], [omega_y], [omega_z]])

    # Substitude Joint Angles Jacobian Matrix
    J_Applied = J.subs({th1: q1, th2: q2, th3: q3, th4: q4, th5: q5, th6: q6})

    # J Transpose
    J_T = J_Applied.transpose()

    # Calculate q_dot
    q_dot = (J_Applied.inv() * X_dot).evalf()

    # Calculate Instantaneous Gravity and Torques
    G_inst = Gravity.subs(
        {th1: q1, th2: q2, th3: q3, th4: q4, th5: q5, th6: q6})

    Taus = G_inst - J_T*F
    Torque1.append(Taus[0])
    Torque2.append(Taus[1])
    Torque3.append(Taus[2])
    Torque4.append(Taus[3])
    Torque5.append(Taus[4])
    Torque6.append(Taus[5])

    # Extracting Elements of q_dot
    q1_dot = q_dot[0]
    q2_dot = q_dot[1]
    q3_dot = q_dot[2]
    q4_dot = q_dot[3]
    q5_dot = q_dot[4]
    q6_dot = q_dot[5]

    # Numerical Integration of Each Joint Angle
    q1 = q1 + q1_dot * dt
    q2 = q2 + q2_dot * dt
    q3 = q3 + q3_dot * dt
    q4 = q4 + q4_dot * dt
    q5 = q5 + q5_dot * dt
    q6 = q6 + q6_dot * dt

## -------------------Plotting---------------------##

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
fig.suptitle('Torques of Joints 1-6')

ax1.plot(time, Torque1)
ax1.set_title("Joint 1")
ax1.set(xlabel='Time (s)')
ax1.set(ylabel='Torque (N-m)')

ax2.plot(time, Torque2)
ax2.set_title("Joint 2")
ax2.set(xlabel='Time (s)')
ax2.set(ylabel='Torque (N-m)')

ax3.plot(time, Torque3)
ax3.set_title("Joint 3")
ax3.set(xlabel='Time (s)')
ax3.set(ylabel='Torque (N-m)')

fig.tight_layout()

ax4.plot(time, Torque4)
ax4.set_title("Joint 4")
ax4.set(xlabel='Time (s)')
ax4.set(ylabel='Torque (N-m)')

ax5.plot(time, Torque5)
ax5.set_title("Joint 5")
ax5.set(xlabel='Time (s)')
ax5.set(ylabel='Torque (N-m)')

ax6.plot(time, Torque6)
ax6.set_title("Joint 6")
ax6.set(xlabel='Time (s)')
ax6.set(ylabel='Torque (N-m)')

fig.tight_layout()


ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax6.grid()
plt.show()
plt.show()
