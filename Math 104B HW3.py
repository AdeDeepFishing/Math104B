'''
Problem 1

(a)

$P_n(x)=f[x_0]+f[x_0,x_1](x-x_0)+f[x_0,x_1,x_2](x-x_0)(x-x_1)+...+f[x_0,...,x_n](x-x_0)...(x-x_{n-1})$   [1]

$      =l_0(x)f_0+...+l_n(x)f_n$  [2]

Now we find the coefficient of $x^n:$

From [1] we know that coefficient of $x^n$ us $f[x_0,...,x_n]$,

so $P_n(x)=l_0(x)f_0+...+l_n(x)f_n=\frac{\prod (x-x_k)}{\prod (x_0-x_k)}f_0+...+\frac{\prod (x-x_k)}{\prod (x_n-x_k)}f_n=\sum_{i=0}^{n}\frac{\prod (x-x_k)}{\prod (x_i-x_k)}f_i $.

And then we have $c_nx^n=\sum_{i=0}^{n}\frac{x^n}{\prod_{k \neq j}(x_j-x_k)}f_i$, so $c_n=\sum_{i=0}^{n}\frac{f_i}{\prod(x_i-x_j)}$.

Therefore, we can have $f[x_0,x_1,...,x_n]=\sum_{j=0}^{n}\frac{f(x_j)}{\prod_{k=0,k \neq j}^{n}(x_j-x_k)}$.


(b)

First, interchange the ith and jth terms, so now we change the position of $\frac{f(x_j)}{\prod_{k=0,k \neq j}^{n}(x_j-x_k)}$ and $\frac{f(x_i)}{\prod_{k=0,k \neq i}^{n}(x_i-x_k)}$. However, the sum of the whole thing would still be the same. 

Thus, $f[x_0,...,x_n]$ stay the same as before. Therefore, the corresponding divided difference unchanged.

'''

#problem 2

import matplotlib.pyplot as plt
import numpy as np

def newton_coeff(x, c, n):
  #return Newton's divided difference
  for i in range(1, n):
    for j in range(n-i):
      c[j][i]=((c[j][i-1]-c[j+1][i-1])/(x[j]-x[i+j]))
  return c

def applyNewtonDD(val, x, c, n):
  #applyr Newton's divided difference's sum
  sum = c[0][0]
  for i in range(1, n):
    temp = 1
    for j in range(i):
      temp = temp * (val - x[j])
    sum = sum + (temp * c[0][i])
  return sum

n = 11
j=np.arange(n)
x=-1+j*2/10
j=np.arange(101)
x_n=-1+j*2/100
y = [[0 for i in range(100)]
        for j in range(100)]

for i in range(n):
  #function for e^(-x^2)
  y[i][0]=np.exp(-x[i]**2)

y=newton_coeff(x, y, n)
y1=[]
for i in range(len(x_n)):
  y1.append(applyNewtonDD(x_n[i], x, y, n))
y2=np.exp(-x_n*x_n)
y1=np.array(y1)

plt.plot(x_n,y2-y1)
plt.show()



'''
Question 3:

First, we know $\bar{x}=f^{-1}(0)$, and $P_1(0)=x_0+f^{-1}[y_0,y_1](0-y_0)=x_0-y_0*f^{-1}[y_0,y_1]$.

Then, applying lagrange's Interpolation, we have $\bar{x}=\frac{y-y_1}{y_0-y_1}x_0+\frac{y-y_0}{y_1-y_0}x_1$.

Now we plug in $y=0$, $y_0=f(0.5)=−0.106530659712633$, $y_1=f(0.6)=0.051188363905973$, $x_0=0.5$, $x_1=0.6$, we have $\bar{x}=\frac{y-y_1}{y_0-y_1}x_0+\frac{y-y_0}{y_1-y_0}x_1=\frac{0-0.051188363905973}{−0.106530659712633-0.051188363905973}0.5+\frac{0-−0.106530659712633}{0.051188363905973-−0.106530659712633}0.6$=$0.5675445848373015508226692277430090326060931430933042561620701094$.
'''

'''
Question 4:

By the Hermite interpolation polynomial, we have $H_3(x)=(1+2\frac{x-x_0}{x_1-x_0}))(\frac{x_1-x}{x_1-x_0})^2*f(x_0)+(x-x_0)(\frac{x_1-x}{x_1-x_0})^2*f'(x_0)+(1+2\frac{x_1-x}{x_1-x_0})(\frac{x_0-x}{x_0-x_1})^2*f(x_1)+(x-x_1)(\frac{x_0-x}{x_0-x_1})^2*f'(x_1)=(1+2\frac{x-0}{1-0})(\frac{1-x}{1-0})^2*0+(x-0)(\frac{1-x}{1-0})^2*0+(1+2\frac{1-x}{1-0})(\frac{0-x}{0-1})^2*2+(x-1)(\frac{0-x}{0-1})^2*3=(1+2(1-x))*2x^2+(x-1)*3x^2=2x^2+4x^2-4x^3+3x^3-3x^2=-x^3+3x^2$.
'''

#problem 5&6

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

def jacobiIterative(A,b,N=25,x=None):
    #Using Jacobi iterative to solve the equation Ax=b
    #https://en.wikipedia.org/wiki/Jacobi_method
    if x is None:
        x = np.zeros(len(A[0]))    # prepare x with all zeros first                                                                                                                                                                
    D = np.diag(A)
    R = A - np.diagflat(D)                                                                                                                                                                     
    for i in range(N):
        x = (b - np.dot(R,x))/D
    return x
  
  def cubicSpline(x, y, tol = 1e-100):
  #return the array lists of b,c,d
  x = np.array(x)
  y = np.array(y)
  if np.any(np.diff(x) < 0):
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
  size = len(x)
  deltaX = np.diff(x)
  deltaY = np.diff(y)
  A = np.zeros(shape = (size,size))  #matrix A
  b = np.zeros(shape=(size,1))  #matrix b
  A[0,0] = 1
  A[-1,-1] = 1

  for i in range(1,size-1):
    #set up the matrixs here
    A[i, i-1] = deltaX[i-1]
    A[i, i+1] = deltaX[i]
    A[i,i] = 2*(deltaX[i-1]+deltaX[i])
    b[i,0] = 3*(deltaY[i]/deltaX[i] - deltaY[i-1]/deltaX[i-1])

  c=np.linalg.solve(A,b)
  d = np.zeros(shape = (size-1,1))
  b = np.zeros(shape = (size-1,1))
  for i in range(0,len(d)):
    d[i] = (c[i+1] - c[i]) / (3*deltaX[i])
    b[i] = (deltaY[i]/deltaX[i]) - (deltaX[i]/3)*(2*c[i] + c[i+1])

  #squeeze() function: remove single-dimensional entries from the shape of an array.
  return b.squeeze(), c.squeeze(), d.squeeze()

#testing
x=[0,0.2,0.4,0.6,0.8]
y=[0.7,0.2,0.1,0.6,1.2]
(b,c,d)=cubicSpline(x,y,1e-100)
print("Test: the coefficients are:")
print('b:',b)
print('c:',c)
print('d:',d)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ax = plt.figure().gca(projection='3d')
#data from #6:
t=[0,.618,.935,1.255,1.636,1.905,2.317,2.827,3.33]
x=[1.5,.9,.6,.35,.2,.1,.5,1,1.5]
y=[.75,.9,1,.8,.45,.2,.1,.2,.25]

(b,c,d)=cubicSpline(x,y,1e-100)
(bx,cx,dx)=cubicSpline(x,t,1e-100)
(by,cy,dy)=cubicSpline(y,t,1e-100)

# #actually i have no idea how to cooporate with #5 here
# ax.plot(b,c,d)
# ax.plot(bx,cx,dx)
# ax.plot(by,cy,dy)

ax.plot(x,y,t)

plt.show()







