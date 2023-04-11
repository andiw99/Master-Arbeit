import sympy as sp

def main():
    t = sp.symbols("t")
    x, v, xv, = sp.symbols("x v xv", cls=sp.Function)
    k, m, eta, D = sp.symbols("k m \eta D", positive=True)

    eqs = [
        sp.Eq(x(t).diff(t), v(t)),
        sp.Eq(v(t).diff(t), -x(t))
    ]

    eqs = [
        sp.Eq(x(t).diff(t), 2 * xv(t)),
        sp.Eq(xv.diff(t), v(t) - k/m * x(t) - eta/m * xv(t)),
        sp.Eq(v.diff(t), 2*D - 2 * eta/m * v(t) - 2* k/m * xv(t))
    ]
    """
    t0 = 0
    x0 = 0
    v0 = 0
    xv0 = 0

    v0 = sp.symbols("v0")
    xv0 = sp.symbols("xv0")


    ics = {
        x(t0): x0,
        v(t0): v0,
        xv(t0): xv0,
    }
    """
    sol=sp.dsolve(eqs,  [x(t), v(t), xv(t)])
    print(sol)
    print(sp.latex(sol))

if __name__ == "__main__":
    main()