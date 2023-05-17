#Porblem 1
import numpy as np

def tridiagonal_solver(a, b, c, d):
    """
    Solves the tridiagonal system of equations Ax = d by using the Thomas algorithm.
    I know we should have Ax=b, but we already use b for other thing, so d here.
    for the input:
      a: array containing the subdiagonal entries of A (length n-1)
      b: array containing the diagonal entries of A (length n)
      c: array containing the superdiagonal entries of A (length n-1)
      d: array containing the right-hand side of the system (length n)
    return: array containing the solution to the system (length n)
    """
    n=len(b)
    c_star=np.zeros(n-1)
    d_star=np.zeros(n)

    # Forward elimination
    c_star[0]=c[0] / b[0]
    d_star[0]=d[0] / b[0]
    for i in range(1, n-1):
        c_star[i]=c[i] / (b[i] - a[i-1]*c_star[i-1])
    for i in range(1, n):
        d_star[i]=(d[i] - a[i-1]*d_star[i-1]) / (b[i] - a[i-1]*c_star[i-1])

    # Backward substitution
    x=np.zeros(n)
    x[-1]=d_star[-1]
    for i in range(n-2, -1, -1):
        x[i]=d_star[i] - c_star[i]*x[i+1]

    return x

  
  #testing

#generate a random tridiagonal matrix A and a corresponding vector x
n=6  # size of matrix, we could change
a=np.random.rand(n-1)
b=np.random.rand(n)
c=np.random.rand(n-1)
x=np.arange(1, n+1)
#x=np.random.rand(n)
A=np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1)
d=np.dot(A, x) # Ax=d

print("The orginal x is")
print(x)
print(" ")

#temp
print("a is:")
print(a)
print("b is:")
print(b)
print("c is:")
print(c)
print(" ")
#tempEnd


print("The matrix A is")
print(A)
print("We have d, which is")
print(d)
print(" ")

x_new=tridiagonal_solver(a, b, c, d)
print("We now have solution from the solver: ")
print(x_new)
print(" ")
#check if the computed solution vector x_new is close enough to the original solution vector x
print("Is the answer the same as our orginal x?")
print(np.allclose(x, x_new)) 



#problem 2
import numpy as np
import matplotlib.pyplot as plt
print("(a)")

N=50
h=1/N

a=[]
b=[]
c=[]
d=[]

#x_j=np.linspace(0,1,51)

for i in range(1,N+1):
  x_i=(1/(N+1))*i
  b.append((2/(h*h))+np.pi**2)
  d.append(2*(np.pi**2)*np.sin(np.pi*(x_i)))

for i in range(1,N):
  a.append((-1/(h*h)))
  c.append((-1/(h*h)))

# print(a)
# print(b)
# print(c)
# print(d)

# A=np.diag(a, -1) + np.diag(b, 0) + np.diag(c, 1)
# print(A)

x_new=tridiagonal_solver(a, b, c, d)
x_newnew=[]
x_newnew.append(0)
for i in range(N):
  x_newnew.append(x_new[i])
x_newnew.append(0)
print("we have x = ")
print(x_newnew)



plt.plot(np.linspace(0,1,N+2),np.array(x_newnew))
plt.show()

print("(b)")
u=[]

for x in np.linspace(0,1,N+2):
  u.append(np.sin(np.pi*x))

plt.plot(np.linspace(0,1,N+2),u,label="Actual")
plt.plot(np.linspace(0,1,N+2),np.array(x_newnew),label="Approximation")
plt.legend()
plt.show()

print("(c)")

totalNotRootYet=0
for i in range(N+2):
  temp=(x_newnew[i]-u[i])**2
  totalNotRootYet+=temp
total=np.sqrt(totalNotRootYet)
print("The error of our approximation in the 2-norm when N=50 is: ")
print(total)


print("Now we take N=100")
N=100
h=1/N

a=[]
b=[]
c=[]
d=[]

#x_j=np.linspace(0,1,51)

for i in range(1,N+1):
  x_i=(1/(N+1))*i
  b.append((2/(h*h))+np.pi**2)
  d.append(2*(np.pi**2)*np.sin(np.pi*(x_i)))

for i in range(1,N):
  a.append((-1/(h*h)))
  c.append((-1/(h*h)))

x_new=tridiagonal_solver(a, b, c, d)
x_newnew1=[]
x_newnew1.append(0)
for i in range(N):
  x_newnew1.append(x_new[i])
x_newnew1.append(0)
print("we have x = ")
print(x_newnew1)

u1=[]

for x in np.linspace(0,1,N+2):
  u1.append(np.sin(np.pi*x))

plt.plot(np.linspace(0,1,N+2),u1,label="Actual")
plt.plot(np.linspace(0,1,N+2),np.array(x_newnew1),label="Approximation")
plt.legend()
plt.show()


totalNotRootYet1=0
for i in range(N+2):
  temp=(x_newnew1[i]-u1[i])**2
  totalNotRootYet1+=temp
total1=np.sqrt(totalNotRootYet1)
print("The error of our approximation in the 2-norm when N=100 is: ")
print(total1)

print("I expected the error would be 1/4. However we can see:")
print(total1/total)
print("Therefore, there are other factors that affect our result so that it does not show what we wanted 🙂")

'''
(d)

If we don't know the exact solution, we can test how different the numbers are. To be specifically, I mean we can have N=20, N=40,N=1000, and N=2000. We can tell that the difference from 20 to 40 is much bigger than the one from 1000 to 2000, which means we can tell the result by looking at how fast the difference they are changing.
'''
