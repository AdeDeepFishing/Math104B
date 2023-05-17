'''
Question 1:

(a) & (b)
'''

import numpy as np
import sys

#get the input of unknown numbers x
n=int(input("Enter how many of unknown (x) we need to solve: "))

a=np.zeros((n,n+1))

#make an array with size n and make them zero first to store solution vector later
x=np.zeros(n)

#get the matrix coefficients
print("\nEnter the Matrix Coefficients details:")
for i in range(n):
    for j in range(n+1):
       #get each of the coefficient, left part is A, and the rightest would be b, for Ax=b
        a[i][j]=float(input( "a["+str(i+1)+"]["+ str(j+1)+"]="))

#Gauss Elimination
for i in range(n):
    if a[i][i]==0.0:
        sys.exit('Oh no, we found there is something divide by zero 😵!')
    for j in range(i+1, n):
        temp=a[j][i]/a[i][i]
        for k in range(n+1):
            a[j][k]=a[j][k]-temp*a[i][k]

#get the value of x
x[n-1]=a[n-1][n]/a[n-1][n-1]
for i in range(n-2,-1,-1):
    x[i]=a[i][n]
    for j in range(i+1,n):
        x[i]=x[i]-a[i][j]*x[j]
    x[i]=x[i]/a[i][i]

print('\nTherefore, the solution is: ')
for i in range(n):
    print('X%d = %0.2f' %(i,x[i]), end = '\t')

#(c)
import numpy as np
import sys

#get the input of unknown numbers x
n=int(input("Enter how many of unknown (x) we need to solve: "))

a=np.zeros((n,n+1))

#make an array with size n and make them zero first to store solution vector later
x=np.zeros(n)

#get the matrix coefficients
print("\nEnter the Matrix Coefficients details:")
for i in range(n):
    for j in range(n+1):
       #get each of the coefficient, left part is A, and the rightest would be b, for Ax=b
        a[i][j]=float(input( "a["+str(i+1)+"]["+ str(j+1)+"]="))

#Gauss Elimination
for i in range(n):
    if a[i][i]==0.0:
        sys.exit('Oh no, we found there is something divide by zero 😵!')
    for j in range(i+1, n):
        temp=a[j][i]/a[i][i]
        for k in range(n+1):
            a[j][k]=a[j][k]-temp*a[i][k]

#get the value of x
x[n-1]=a[n-1][n]/a[n-1][n-1]
for i in range(n-2,-1,-1):
    x[i]=a[i][n]
    for j in range(i+1,n):
        x[i]=x[i]-a[i][j]*x[j]
    x[i]=x[i]/a[i][i]

print('\nTherefore, the solution is: ')
for i in range(n):
    print('X%d' %(i,x[i]), end = '\t')
