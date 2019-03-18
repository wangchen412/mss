from subprocess import Popen, PIPE
from scipy.optimize import minimize
from scipy.optimize import basinhopping


p = Popen(['./mismatch'], shell=True, stdout=PIPE, stdin=PIPE)
def func(x):
    p.stdin.write(bytes("".join([str(i) for i in x])[:-1] + '\n', 'UTF-8'))
    p.stdin.flush()
    return float(p.stdout.readline().strip())

print(basinhopping(func, [1, 1, 1, 1],
                   minimizer_kwargs={"method":"Nelder-Mead", "tol":1e-4},
                   disp=True, niter=20))
