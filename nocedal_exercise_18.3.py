import scipy.optimize as opt
import numpy as np
import sympy as sy

# firstly, use package.
x0 = np.array([-1.71, 1.59, 1.82, -0.763, -0.763])


def obj_fun(x):
    return np.exp(x[0] * x[1] * x[2] * x[3] * x[4]) - 0.5 * (x[0] ** 3 + x[1] ** 3 + 1) ** 2


cons = [{'type': 'eq', 'fun': lambda x: x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 + x[4] ** 2 - 10},
        {'type': 'eq', 'fun': lambda x: x[1] * x[2] - 5 * x[3] * x[4]},
        {'type': 'eq', 'fun': lambda x: x[0] ** 3 + x[1] ** 3 + 1}
        ]

res = opt.minimize(obj_fun, x0, method='SLSQP', constraints=cons, tol=1e-6)

####################################################################################

# secondly, write Algorithm 18.1 from Nocedal

# symbols init
x = sy.symbols('x[0], x[1], x[2], x[3], x[4]')
l = sy.symbols('l[0] l[1] l[2]')
m = len(l)
n = len(x)
convergence_tol = 1e-6

# functions are init
obj_fuc_expr = sy.exp(x[0] * x[1] * x[2] * x[3] * x[4]) - 0.5 * (x[0] ** 3 + x[1] ** 3 + 1) ** 2
con = (x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 + x[4] ** 2 - 10, x[1] * x[2] - 5 * x[3] * x[4],
       x[0] ** 3 + x[1] ** 3 + 1)
lagranian = obj_fuc_expr + np.dot(l, con)

# derivatives of functions
lagrangian_xx = sy.Matrix([[sy.diff(lagranian, x[i], x[j]) for i in range(n)] for j in range(n)])
a = sy.Matrix([[sy.diff(con[i], x[j]) for j in range(n)] for i in range(m)])

# coefficient matrix construction
coefficient = sy.zeros(m + n, m + n)
coefficient[:n, :n] = lagrangian_xx
coefficient[:n, n:] = - a.T
coefficient[n:, :n] = a

# rhs of the linear eq
b = sy.Matrix([sy.Matrix([sy.diff(-obj_fuc_expr, x[i]) for i in range(n)]), -sy.Matrix(con)])

F = sy.Matrix([sy.diff(lagranian, x[i]) for i in range(n)] + list(con))
Fn = sy.lambdify((x, l), F, "numpy")
# here error is defined as L2 norm of F function
error = lambda x_temp, l_temp: np.dot(Fn(x_temp, l_temp).flatten(), Fn(x_temp, l_temp)[:, 0].flatten())


def numerical_matrix(A, b, x_k, l_k):
    """
    plugging in the numbers at the point (x_k, l_k) to A matrix and b column vector
    Args:
        A (sy.Matrix): coefficient matrix
        b (sy.Matrix): rhs column
        x_k (np.ndarray): point x_k
        l_k (np.ndarray): estimate of l

    Returns:
        p (np.ndarray): sympy matrix to denote direction

    """
    A_k = sy.lambdify((x, l), A, "numpy")
    b_k = sy.lambdify((x,), b, "numpy")

    return A_k(x_k, l_k), b_k(x_k)


# initialisation of the SQP method
a_f = sy.lambdify((x,), a)
a_k = a_f(x0)  # type: np.ndarray
df = sy.lambdify((x,), b[:n, 0])
delf = df(x0)
l0 = np.matmul(np.matmul(np.linalg.inv(np.matmul(a_k, a_k.T)), a_k), delf).flatten()
x_k, l_k = x0, l0

while error(x_k, l_k) > convergence_tol:
    A_k, b_k = numerical_matrix(coefficient, b, x_k, l_k)
    print("\n\ndeterminant of coefficient matrix", np.linalg.det(A_k))
    print("x_k = ", x_k)
    print("l_k = ", l_k)
    p = np.linalg.solve(A_k, b_k).flatten()
    p_k = p[:n]
    l_k1 = p[n:]
    x_k = x_k + p_k
    l_k = l_k1
    print("p = ", p)

f = obj_fuc_expr.subs(zip(x, x_k))

print("the solution", x_k, "\nwith obj_func_val", f)

# the second method does not converges. sensitive to the choice of l0.