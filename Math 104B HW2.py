import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from collections.abc import Callable
from typing import List,Tuple
import math

'''
Problem 1:

(a)

We first find the Lagrange basis polynomials for the data set $x = [0,1,3]$ and $y = [1,1,-5]$.

$l_1(x)=\frac{(x-x_2)(x-x_3)}{(x_1-x_2)(x_1-x_3)}=\frac{(x-1)(x-3)}{(0-1)(0-3)}= \frac{1}{3}x^{2}-\frac{4}{3}x+1$

$l_2(x)=\frac{(x-x_1)(x-x_3)}{(x_2-x_1)(x_2-x_3)}=\frac{(x-0)(x-3)}{(1-0)(1-3)}= -\frac{1}{2}x^{2}+\frac{3}{2}x$

$l_3(x)=\frac{(x-x_1)(x-x_2)}{(x_3-x_1)(x_3-x_2)}=\frac{(x-0)(x-1)}{(3-0)(3-1)}= \frac{1}{6}x^{2}-\frac{1}{6}x$

Now we can calculate $P_2(x)$:

$P_2(x)=l_1(x)*f_1+l_2(x)*f_2+l_3(x)*f_3=(\frac{1}{3}x^{2}-\frac{4}{3}x+1)*1+(-\frac{1}{2}x^{2}+\frac{3}{2}x)*1+(\frac{1}{6}x^{2}-\frac{1}{6}x)*(-5)=-x^2+x+1$

Therefore, the Lagrangian form of the interpolating polynomial $P_2(x)=-x^2+x+1$


(b)

From (a), we can have $f(2)=P(2)=-(2)^2+(2)+1=-1$
'''

#Problem 2

def BarycentricFormula(x):
    #x means the array of nodes x_0, x_1, ... that input
    #return lambda which is the barycentric weight

    N=x.size
    lmbda=np.zeros(N)
    for j in range(N):
        #construct lambda of j
        product=1
        for k in range(N):
            if(k!=j):
                product*=(x[j]-x[k])
        lmbda[j] = 1/product
    return lmbda
    
x=np.array([1,2,3])
BarycentricFormula(x)

def Poly(x,xData,fData):
    #x is the that we plug in, xData & fData are the points
    
    lmbda=BarycentricFormula(xData)
    N=xData.size
    
    #define the top of Pn(x)
    #summand=np.zeros(N)
    top=0
    #define the bottom of Pn(x)
    #summand1=np.zeros(N)
    down=0
    
    for j in range(N):
        top+=(lmbda[j]*fData[j])/(x-xData[j])
        down+=lmbda[j]/(x-xData[j])
        # summand[j]=(lmbda[j]*fData[j])/(x-xData[j])
        # top+=summand[j]
        # summand1[j]=lmbda[j]/(x-xData[j])
        # down+=summand1[j]
    
    return top/down

xData=np.array([0,0.25,0.5,0.75,1.25,1.5])
fData=np.array([0,0.7071,1,0.7071,-0.7071,-1])
print(Poly(2,xData,fData))


#problem 3
for n in [4,8,12]:
    x_jData=np.zeros(n+1)
    fx_jData=np.zeros(n+1)

    def f1(x):
        return 1/(1+x**2)

    for j in range(n+1):
        x_jData[j] = -5+j*(10/n)

    for j in range(n+1):
        fx_jData[j] = f1(x_jData[j])

    domain=np.linspace(-5,5,5000)
    
    plt.plot(domain,Poly(domain,x_jData,fx_jData))

plt.plot(np.linspace(-5,5,5000),f1(np.linspace(-5,5,5000))) #red line is f(x)
plt.show()


for n in [4,8,12,100]:
    x_jData=np.zeros(n+1)
    fx_jData=np.zeros(n+1)

    def f1(x):
        return 1/(1+x**2)

    for j in range(n+1):
        x_jData[j] = 5*math.cos(j*math.pi/n)

    for j in range(n+1):
        fx_jData[j] = f1(x_jData[j])

    domain=np.linspace(-5,5,5000)
    plt.plot(domain,Poly(domain,x_jData,fx_jData))
plt.plot(np.linspace(-5,5,5000),f1(np.linspace(-5,5,5000))) #purple line is f(x)
plt.show()


for n in [4,8,12]:
    x_jData=np.zeros(n+1)
    fx_jData=np.zeros(n+1)

    def f2(x):
        return np.exp((-x**(2))/5)

    for j in range(n+1):
        x_jData[j] = -5+j*(10/n)

    for j in range(n+1):
        fx_jData[j] = f2(x_jData[j])

    domain=np.linspace(-5,5,5000)
    plt.plot(domain,Poly(domain,x_jData,fx_jData))
plt.plot(np.linspace(-5,5,5000),f2(np.linspace(-5,5,5000))) #red line is f(x)
plt.show()


'''
We can see in all cases, when n is larger, it would be more similar/close to the original function.
'''

#probelm 4

'''
Proof: 

Let $f \in C^2 [x_0,x_1]$ and P1 its interpolation linear polynomial at $x_0$ and $x_1$. Let $|f''(x)|<= M_2$ for all $x \in [x_0,x_1]$. Also, let $||f-P_1||_\infty=max|f(x)-P_1|$, where $x \in [x_0,x_1]$. 

Claim:$f(x) - P_1(x) = \frac{f^{(2)}(\xi)}{2} \prod_{i=0}^n (x-x_i)$.
Proof:We first have the error term $R_1(x) = f(x) - P_1(x)$, and then an auxiliary function would be $Y(t) = R_1(t) - \frac{R_1(x)}{W(x)} W(t) $, with $W(t) = \prod_{i=0}^1 (t-x_i) $. Since $x_i$ are roots of $R_1(t)$ and $W(t)$, we have $Y(x)=Y(x_i)=0$, so we know $Y$ has at least $3$ roots while $n=1$. By Rolle's theorem, we know $Y'(t)$ has at least $n+1$ roots, then $Y''(t)$ has at least one root $ξ$ such that $ξ$ is in the interval $I$. Then, $Y''(t) = R_1''(t) - \frac{R_1(x)}{W(x)} \ 2! $. Also, since $p_n(x)$ is a polynomial, we know $R_1''(t) = f''(t)$. Therefore, $Y''(t) = f''(t) - \frac{R_1(x)}{W(x)} \ 2!$ Since $ξ$ is the root of $Y''(t)$, we have $Y''(ξ)=f''(ξ) - \frac{R_1(x)}{W(x)} \ 2! = 0$. Therefore, $R_1(x) = f(x) - p_1(x) = \frac{f''(ξ)}{2!} \prod_{i=0}^1 (x-x_i)$.

Now we back to our original proof, we know $|x-x_0||x-x_1|<\frac{|x_1-x_0|}{2}*\frac{|x_0-x_1|}{2}$. Then, $|R_1(x)|<= \frac{h^2}{4*(1+1)}max|f''(ξ)|$, with $ξ \in [x_0,x_1]$. Therefore, we have $||f-P_1||_\infty=\frac{1}{8}(x_1-x_0)^2M2$

'''

