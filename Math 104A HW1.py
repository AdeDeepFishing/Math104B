import numpy as np

#Question 2 here:
# have a function of Trapezoidal Rule that take four values: a(lower),b(upper),N(number of intevals bewteen a & b), f(the function)
def trapezoidalRule(a,b,N,f):
    
    #calculate h, the step size
    h=(b-a)/N
    
    #calculate the sum:
    
    #sum of the first value and the last value
    integration=(1/2)*(f(b)+f(a))
    
    #sum of the middle values except the fisrt and the last value
    for i in range(1,int(N)):
        #the next temporary point
        temp=a+i*h
        #add the value to our existing integration
        integration+=f(temp)
    
    #final sum
    integration*=h
    
    #return the value
    return integration

  
  #Question 3 here:

#As the qusetion requires, it is the function we gonna test: f(x)=xe^(x^2)
def f(x):
    return x*np.exp(x**2)

#calculate the ingegral of the function above: (e-1)/2
integralf = (np.exp(1)-1)/2

# for h = 1/10
print("Applying the function of Trapezoidal Rule, while h=1/10, we have",trapezoidalRule(0,1,10,f))
# for h = 1/20
print("Applying the function of Trapezoidal Rule, while h=1/20, we have",trapezoidalRule(0,1,20,f))
# for h = 1/40
print("Applying the function of Trapezoidal Rule, while h=1/40, we have",trapezoidalRule(0,1,30,f))

#calculate the error
error0=integralf-trapezoidalRule(0,1,10,f)
print("We have the error while h=1/10 : ",error0)
error1=integralf-trapezoidalRule(0,1,20,f)
print("We have the error while h=1/20 : ",error1)
error2=integralf-trapezoidalRule(0,1,40,f)
print("We have the error while h=1/40 : ",error2)

#check if T_h has a convergent trend at the expected quadratic rate
if (round(error0/error1)==4):
    print("Since the error of h=1/10 is four times of the error of h=1/20, T_h has a convergent trend at the expected quadratic rate from h=1/10 to h=1/20")
else:
    print("Unfortunately, T_h does not have a convergent trend at the expected quadratic rate from h=1/10 to h=1/20")

if (round(error1/error2)==4):
    print("Since the error of h=1/20 is four times of the error of h=1/40, T_h has a convergent trend at the expected quadratic rate from h=1/20 to h=1/40")
else:
    print("Unfortunately, T_h does not have a convergent trend at the expected quadratic rate from h=1/20 to h=1/40")



#Question 4 here:

#The integral of $\int_{a}^b f(x)dx$ can be approximated by:

def g(x):
    return np.exp(-x**2)

#have two empty lists first, gonna fill in with numbers later
tempArraySize=14
T=np.zeros(tempArraySize)
q=np.zeros(tempArraySize-2)

#calculating all the T in the TList
for i in range(0,tempArraySize):
    N=2**(i)
    T[i]=trapezoidalRule(0,1,N,g)
    
#calculating all the q in the qList
for i in range(0,tempArraySize-2):
    q[i]=(T[i+1]-T[i])/(T[i+2]-T[i+1])

#got the qList
print(q)
print()

#print the one that cloest to 4
print("q(h) is approximately equal to 4: ",q[10])

'''
We can see when $ N=2^{10} $, q is very close to 4, so we can have $ N=2^{10} $ and we will proceed with it in the following.

Try approximate $I-T_h$ for $ N=2^{10} $. We can use $ N=2^{20} $ or anything much greater than N

We can now see the approximation for the error, $I-T_h$ for $N=2^{10}$ is
'''

#calculate the I that we mentioned above
I_true=trapezoidalRule(0,1,2**(20),g)
#calculate the error approximation as required
errorApproximation=I_true-T[10]
print(errorApproximation)

'''
Now we try to apply the error approximation in (b):

$S_h=T_h+\frac{4}{3}(T_{\frac{h}{2}}-T_h)=$
'''
S_h=trapezoidalRule(0,1,2**(10),g)+errorApproximation
print(S_h)

'''
Our initial approximation, $T_h$, can be written in the form $ch^{2}+o(h^{p})$. Applying extrapolation approximation, the error becomes $o(h^{p})$ and it is less than the original error. As $ p \geq 2, o(h^{p})$ approches to zero much faster than the error of inital approximation. Thus, we have $S_h[e^{-x^{2}}]$ is more accurate and converges faster to $I[e^{-x^{2}}]$ than $T_h[e^{-x^{2}}]$.
'''


