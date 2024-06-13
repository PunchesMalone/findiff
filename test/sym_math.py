import sympy as sym
from sympy.core.function import Function
from sympy_utils import make_fort

x1,x2,x3,x4,x5,x6 = sym.symbols('x(1), x(2), x(3), x(4), x(5), x(6)')

ndim = 6
x = [x1, x2, x3,x4,x5,x6]
jac = []
hess = []
onebyr32 = -sym.sqrt(x1**2 + x2**2 + x3**2)**-3
acc = [x4,x5,x6,x1*onebyr32,x2*onebyr32,x3*onebyr32]

# tess = []
for xj in x:
    for xi in acc:
    # Jacobian elements have no symmetry
        thisDeriv =  xi.diff(xj)
        jac.append(thisDeriv)
for xk in x:
    for xj in x:
        for xi in acc:
            thisSecondDeriv =  xi.diff(xj,xk)
            hess.append(thisSecondDeriv)
# for k,xk in enumerate(x):
#     for j,xj in enumerate(x):
#         for i, xi in enumerate(x):
#             thisThirdDeriv =  P.diff(xi,xj,xk)
#             tess.append(thisThirdDeriv)
make_fort(acc,name="acc", prec="wp")
make_fort(jac,name="jac", prec="wp",shape=[6,6])
make_fort(hess,name="hes", prec="wp",shape=[6,6,6])
# make_fort(tess,name="tes", prec="wp",shape=[3,3,3])
