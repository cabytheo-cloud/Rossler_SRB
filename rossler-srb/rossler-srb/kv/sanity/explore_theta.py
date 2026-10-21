# Non-rigorous exploration (double precision, SciPy):
#  * locate the fixed point of P_c on {x=0, x'>0} at cmid,
#  * check that the polar angular speed theta' = (x^2 + a x y + y^2 + y z)/(x^2+y^2)
#    stays positive (and how small it gets) along
#       - the periodic orbit through p_c,
#       - the 15 iterates of the initial point X0(s_m) = p + s_m e_u,
#       - a cloud of points of the local boxes used in the proof.
#  This only guides the rigorous computation (which re-checks theta' > 0 by
#  enclosing the denominator of the reparametrised field).
import numpy as np
from scipy.integrate import solve_ivp

a = b = 0.2
cmid = 4.50698934631197
E = np.array([[0.9999996961880596, -0.8077389541129416],
              [7.795022697815469e-4, -0.5895403141503822]])
eu, es = E[:, 0], E[:, 1]
p0 = np.array([-6.938698731191516, 0.03606348866252088])
sm = 5.3682155586e-4


def rhs(t, X, c):
    x, y, z = X
    return [-(y + z), x + a * y, b + z * (x - c)]


def thetadot(X):
    x, y, z = X
    return (x * x + a * x * y + y * y + y * z) / (x * x + y * y)


def ev(t, X, c):
    return X[0]
ev.direction = 1
ev.terminal = False


def P(yz, c, dense=False):
    X0 = [0.0, yz[0], yz[1]]
    sol = solve_ivp(rhs, [0, 20], X0, args=(c,), events=ev, rtol=1e-12, atol=1e-14,
                    dense_output=True, max_step=0.01)
    # first positive crossing with t > small
    te = [t for t in sol.t_events[0] if t > 1e-6][0]
    Xe = sol.sol(te)
    ts = np.linspace(0, te, 4000)
    td = np.array([thetadot(sol.sol(t)) for t in ts])
    return np.array([Xe[1], Xe[2]]), td.min(), te


q, m, T = P(p0, cmid)
print("P(p0)-p0 =", q - p0, " min theta' on periodic orbit:", m, " return time", T)

X = p0 + sm * eu
mins = []
for k in range(15):
    X, m, T = P(X, cmid)
    mins.append(m)
print("tip after 15 iterates (U,V) =", np.linalg.solve(E, X - p0))
print("min theta' per iterate:", np.round(mins, 4))

# random points of the W^u box p + E([-6e-4,6e-4]x[-1e-5,1e-5])
rng = np.random.default_rng(1)
mm = 10
for _ in range(20):
    u, v = rng.uniform(-6e-4, 6e-4), rng.uniform(-1e-5, 1e-5)
    _, m, _ = P(p0 + E @ np.array([u, v]), cmid)
    mm = min(mm, m)
print("min theta' over W^u box samples:", mm)

# theta as function of time along the periodic orbit: where is the minimum?
sol = solve_ivp(rhs, [0, T], [0.0, p0[0], p0[1]], args=(cmid,), rtol=1e-12, atol=1e-14, max_step=0.01)
td = np.array([thetadot(sol.y[:, i]) for i in range(sol.y.shape[1])])
i = td.argmin()
print("argmin theta' on periodic orbit: X =", sol.y[:, i], " theta =", np.arctan2(sol.y[1, i], sol.y[0, i]))
print("max z on periodic orbit:", sol.y[2].max(), " max x:", sol.y[0].max())
