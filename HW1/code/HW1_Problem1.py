## Brendan Neal
## ENPM662 Homework 1 Question 1.1

##------------------------Importing Libraries-------------------------##
import numpy as np
from matplotlib import pyplot as plt

##------------------------Variable Definitions------------------------##

X_i = 0 #Initial X (m)
Y_i = 0 #Initial Y (m)
Theta_i = 0 #Initial Orientation (rad)

omega = 5 #Rear wheel speed: rad/s
r = 0.5 #Wheel radius (m)
L = 1.5 #Track Length (m)


T = 10.1 #Time period: seconds
dt = 0.1 #Delta T
TimeData = np.arange(0, T, dt) #Time vector (for plotting)


x_data = [] #Array to store computed X positions for plotting.
y_data = [] #Array to store computed Y positions for plotting.
theta_data = [] #Array to store computed Theta positions for plotting.

##-------------------Euler (Numeric) Integration----------------------##

x = X_i #Set current X to initial X
y = Y_i #Set current Y to initial Y
theta = Theta_i #Set current Theta to initial Theta


#Begin Euler (Numeric) Integration
for t in TimeData:

    #Append Data for Plotting
    x_data.append(x)
    y_data.append(y)
    theta_data.append(theta)

    #State Equations Derived by Hand
    thetaDot = (omega*r*np.tan(0.5*np.sin(np.pi * t)))/L
    xDot = omega*r*np.cos(theta + 0.5*np.sin(np.pi * t))
    yDot = omega*r*np.sin(theta + 0.5*np.sin(np.pi * t))

    #Numerical Integration
    theta = theta + thetaDot*dt
    x = x + xDot*dt
    y = y + yDot*dt

##------------------Plotting Results-------------------------##

figure, (ax1, ax2, ax3) = plt.subplots(1,3)

ax1.set_title('X Position vs Time')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('X Position (m)')
ax1.plot(TimeData, x_data, 'b-')

ax2.set_title('Y Position vs Time')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Y Position (m)')
ax2.plot(TimeData, y_data, 'r-')

ax3.set_title('X Position vs Y Position')
ax3.set_xlabel('X Position (m)')
ax3.set_ylabel('Y Position (m)')
ax3.plot(x_data, y_data, 'ko')

plt.show()


