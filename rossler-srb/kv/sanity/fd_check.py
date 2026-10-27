# Non-rigorous sanity check (double precision, SciPy) of the jet layer of the Kv program:
# first and second derivatives of the local map f(u,v) = E^{-1}(P(p + E(u,v)) - p) at (0,0),
# by central finite differences, with P computed in CARTESIAN coordinates by event location
# on {x = 0, x' > 0} (independent of the polar formulation used in the rigorous code).
import numpy as np
from scipy.integrate import solve_ivp

a = b = 0.2
c = 4.50698934631197
E = np.array([[0.9999996961880596, -0.8077389541129416],
              [7.795022697815469e-4, -0.5895403141503822]])
Ei = np.linalg.inv(E)

def rhs(t, X):
    x, y, z = X
    return [-(y + z), x + a * y, b + z * (x - c)]

def ev(t, X):
    return X[0]
ev.direction = 1

def P(yz):
    sol = solve_ivp(rhs, [0, 12], [0.0, yz[0], yz[1]], method='DOP853', events=ev,
                    rtol=1e-13, atol=1e-15)
    te = [k for k, t in enumerate(sol.t_events[0]) if t > 1e-3][0]
    Xe = sol.y_events[0][te]
    return np.array([Xe[1], Xe[2]])

# fixed point by Newton (finite-difference Jacobian)
p = np.array([-6.938698731191516, 0.03606348866252088])
for _ in range(3):
    F = P(p) - p
    h = 1e-7
    J = np.column_stack([(P(p + h * e) - P(p - h * e)) / (2 * h) for e in np.eye(2)]) - np.eye(2)
    p = p - np.linalg.solve(J, F)
print("p =", repr(p[0]), repr(p[1]))

def f(u, v):
    return Ei @ (P(p + E @ np.array([u, v])) - p)

h1 = 1e-6
a1 = (f(h1, 0) - f(-h1, 0)) / (2 * h1)
a2 = (f(0, h1) - f(0, -h1)) / (2 * h1)
print("a  (columns) :", a1, a2)
h2 = 1e-4
f0 = f(0, 0)
h11 = (f(h2, 0) - 2 * f0 + f(-h2, 0)) / h2**2
h22 = (f(0, h2) - 2 * f0 + f(0, -h2)) / h2**2
h12 = (f(h2, h2) - f(h2, -h2) - f(-h2, h2) + f(-h2, -h2)) / (4 * h2**2)
print("h^1 (11,12,22):", h11[0], h12[0], h22[0])
print("h^2 (11,12,22):", h11[1], h12[1], h22[1])
