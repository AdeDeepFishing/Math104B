#Question 3

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import cg

# Define the function to solve the system using the Jacobi method
def jacobi_solver(A, b, x, max_iter=1000, tol=1e-15):
    n = len(b)
    for k in range(max_iter):
        x_new = np.zeros(n)
        for i in range(n):
            x_new[i] = (b[i] - np.dot(A[i,:], x) + A[i,i]*x[i]) / A[i,i]
        if np.linalg.norm(x_new - x) < tol:
            break
        x = x_new
    return x, k+1

# Set up problem parameters
Ns = [50, 100, 200]
x_exact = np.arange(1, Ns[0]+1)

# Solve using Jacobi method and conjugate gradient method
for N in Ns:
    h = 1/(N+1)
    x = np.linspace(h, 1-h, N)
    A = (-1/h**2)*diags([1, -2, 1], [-1, 0, 1], shape=(N, N)).toarray()
    A += np.eye(N)*(np.pi**2)
    b = 2*np.pi**2*np.sin(np.pi*x)
    x_hat, num_iter = jacobi_solver(A, b, np.zeros(N))
    print(f"N = {N}, Jacobi method took {num_iter} iterations")
    x_hat, num_iter = cg(A, b, x0=np.zeros(N), tol=1e-15, maxiter=1000)
    print(f"N = {N}, Conjugate gradient method took {num_iter} iterations")
    rel_res_norm = np.linalg.norm(b - A @ x_hat) / np.linalg.norm(b)
    if rel_res_norm < 1e-15:
        print(f"N = {N}, Conjugate gradient method converged with final relative residual norm = {rel_res_norm}")
    else:
        print(f"N = {N}, Conjugate gradient method did not converge after {num_iter} iterations with final relative residual norm = {rel_res_norm}")

'''
Based on the results, it seems like the conjugate gradient method did not converge after the maximum number of iterations for all values of N. The final relative residual norm for the conjugate gradient method is also shown, which indicates that the residual was still relatively large even after 1000 iterations. In contrast, the Jacobi method took 1000 iterations for all values of N. (btw, 1000, the number I picked was arbitrary.)

It seems like the conjugate gradient method may require more than 1000 iterations to converge for the given problem. This is because the final relative residual norm for the conjugate gradient method was still relatively large even after 1000 iterations. However, it is also possible that adjusting the stopping criteria or using a different preconditioner may improve the convergence of the conjugate gradient method. Overall, the results suggest that the Jacobi method may be more efficient for this particular problem compared to the conjugate gradient method.
'''
