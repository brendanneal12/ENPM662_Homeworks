##------------------------Importing Libraries-------------------------##
import sympy
##--------------------Defining Symbolic Variables---------------------##

#Defining Time
t = sympy.symbols('t') #Time (s)

a = sympy.symbols('a')
b = sympy.symbols('b')

x = sympy.Function('x')(t)
z = sympy.Function('z')(t)

y = a*sympy.sin(x+z) + b

Y_Dot = sympy.diff(y,t)

Y_Dot_Expanded = sympy.expand_trig(Y_Dot)

Y_Dot_Evaled = Y_Dot_Expanded.expand()

X_Coeffs = sympy.collect(Y_Dot_Evaled, x.diff())
Z_Coeffs = sympy.collect(Y_Dot_Evaled, z.diff())

Coeffs_X_Dot = X_Coeffs.coeff(x.diff())

print(X_Coeffs)
print(Z_Coeffs)

print(Coeffs_X_Dot)




