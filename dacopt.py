import scipy as sp
import scipy.linalg as l


def finDiff(f, xstar,eps=1e-5,order=3):
    # Perturb initial conditions for deriv calculation
    bigpertmatrix = sp.array([xstar for _ in xstar])
    pstep = bigpertmatrix + sp.eye(len(xstar))*eps
    mstep = bigpertmatrix - sp.eye(len(xstar))*eps
    pdf = sp.array([f(step) for step in pstep])
    mdf = sp.array([f(step) for step in mstep])

    if order > 3 :
        twopstep = bigpertmatrix + sp.eye(len(xstar))*2*eps
        twomstep = bigpertmatrix - sp.eye(len(xstar))*2*eps
        twopdf = sp.array([f(step) for step in twopstep])
        twomdf = sp.array([f(step) for step in twomstep])
    if order > 5:
        threepstep = bigpertmatrix + sp.eye(len(xstar))*3*eps
        threemstep = bigpertmatrix - sp.eye(len(xstar))*3*eps
        threepdf = sp.array([f(step) for step in threepstep])
        threemdf = sp.array([f(step) for step in threemstep])
    if order > 7:
        fourpstep = bigpertmatrix + sp.eye(len(xstar))*4*eps
        fourmstep = bigpertmatrix - sp.eye(len(xstar))*4*eps
        fourpdf = sp.array([f(step) for step in fourpstep])
        fourmdf = sp.array([f(step) for step in fourmstep])

    if order == 3:
        dfdx = (pdf - mdf) / (2*eps)
    elif order == 5:
        dfdx = (8*pdf - twopdf - 8*mdf + twomdf)/(12*eps)
    elif order == 7:
        dfdx = (45*pdf - 9*twopdf + threepdf - 45*mdf + 9*twomdf -threemdf)/(60*eps)
    elif order == 9:
        dfdx = (3*(fourmdf - fourpdf) + 32*(threepdf-threemdf) + 168*(twomdf-twopdf) + 672*(pdf-mdf))/(840*eps)
    else:
        print('Warning: finite difference order not implemented. Using 3 point stencil')
        dfdx = (pdf - mdf) /(2*eps)
    return dfdx

def cStep(f,xstar,eps=1e-80):
    # Perturb initial conditions
    bigpertmatrix = sp.array([xstar for _ in xstar])
    cstep = bigpertmatrix + sp.eye(len(xstar))*1j*eps
    cdf = sp.array([f(step) for step in cstep])
    dfdx = sp.imag(cdf)/eps

    return dfdx

def STMchain(STMs):
    # takes list of n node-to-node STMs (OR STTs?):
    # in = [STM(0,1), STM(1,2), . . ., STM(n-1,n)]
    # and produces STM(0,n)
    intSTM = sp.eye(STMs[0].shape[0])
    for STM in STMs:
        intSTM = STM@intSTM
    return intSTM

def fcount():
    fcount.count +=1
def gcount():
    gcount.count +=1
def hcount():
    hcount.count +=1

def TRDX(thisg,thish,TR):
    # spectral decompose and sort H
    evals, evecs = l.eig(thish)
    idx = evals.argsort()
    evals = evals[idx]
    V = evecs[:,idx]
    evals = sp.real(evals)
    D = sp.diag(evals)
    gtil = V.T@thisg
    def dx(lam):
        val = sp.sum([-gtil[i]/(l+lam)*V[:,i] for i, l in enumerate(evals)],0)
        return val
    def phi(lam):
        val = 1/sp.sqrt(sp.sum([(gtil[i]/(l+lam))**2 for i, l in enumerate(evals)],0)) - 1/TR
        return val


    if l.norm(dx(0)) < TR and l.det(thish) > 0:
        return dx(0)

    elif l.det(thish) > 0:
        # use secant method to rootsolve for lambda
        lamb = secantsolve(phi, 1, 2, 1e-5)
        delt = dx(lamb)
        return delt

    elif l.det(thish) <= 0 and l.norm(dx(-evals[0]+1e-4)) > TR:

        lamb1 = secantsolve(phi, 1, 2,  1e-8)
        print('!')
        if lamb1 > -min(evals):
            delt = dx(lamb1)
        else:
            lamb2 = secantsolve(phi, -abs(evals[0])+0.000001, abs(evals[0])+0.000002,1e-5)
            delt = dx(lamb2)
        return delt

    else:
        print('Hard case is hard')
        return dx(0)


def TRMin(func,g,h,x0,TR,eps=1e-4,gammai=1.1,gammad=0.75,acc=0.9,fail=0.2,imax=100):
    xs = [x0]
    gs = [g(x0)]
    hs = [h(x0)]
    Js = [func(x0)]
    iflag = 0
    for i in range(imax):
        dx = TRDX(gs[-1],hs[-1],TR)
        dJexp = gs[-1]@dx + 1/2*dx.T@hs[-1]@dx
        Jtest = func(xs[-1]+dx)
        dJreal = Jtest - Js[-1]

        # TR adjustment block
        rhoe = dJreal/dJexp
        if abs(1-rhoe) > 2*fail or dJreal > 0:
            # took a very bad step, reduce trust region and try again
            TR *= gammad**2
            continue

        elif abs(1-rhoe) < acc:
            # took a good step, increase trust region and accept
            TR *= gammai

        elif abs(1-rhoe) < fail:
            # took a mildly bad step, reduce trust region but accept step
            TR *= gammad

        xs.append(xs[-1]+dx)
        gs.append(g(xs[-1]))
        hs.append(h(xs[-1]))
        Js.append(Jtest)
        if l.norm(gs[-1]) < eps:
            break
        if i == imax-1:
            print('Too many iterations')
            iflag = 1
    return [xs, Js, gs, hs, iflag]



def polymin(ffile,x0,eps,stride,tmul,imax,disp=0):
    # Adapted from equivalent MATLAB by Russell, OST class notes
    # ffile is function handle for the 1D function f(x) to minimize
    # x0 is the starting location for x
    # eps is a tolerance for when to stop refining
    # stride is the initial stride length for bracketing the min\
    # tmul is a refining parameter (0<tmul<1) that determines how aggressively to refine
    # imax is the max number of refine iterations
    # dependencies:
    # setBounds() returns 3 points for x and f with min bracketed
    # polyminOnce() takes three points and assoc. f evals and returns the min if it were a
    # quadratic
    iflag=0
    f0=ffile(x0)
    x = sp.array([0, x0, 0]).T
    f = sp.array([0, f0, 0]).T # middle function eval

    tol=1
    for i in range(imax):
        if(i>0):
            stride=stride*tmul
        [x, f, iflag]=setBounds(ffile,x[1],f[1],stride,disp)
        xstar = polyminOnce(x,f)
        fstar = ffile(xstar)
        if disp > 0.01: print('cost =',fstar)
        tol=sp.absolute(max(f)-min(f))/max([1e-8,sp.absolute(min(f))])
        if(f[1]>fstar):
            x[1]=xstar
            f[1]=fstar
        if(i==imax-1):
            print('warning max iteration reached in polymin')
            iflag=1
        if(tol) <= eps:
            break
    return [xstar, fstar, iflag]


def polyminOnce(x,f):
    # Take three ordered pairs as input, fit them to a quadratic, return the min. 
    # y = a*x^2 + b*x + c
    x = sp.array(x)
    f = sp.array(f)
    xmin = (f[0]*x[1]**2 - f[0]*x[2]**2 - f[1]*x[0]**2 +\
            f[1]*x[2]**2 + f[2]*x[0]**2-f[2]*x[1]**2)/\
            (2*(f[0]*x[1] - f[0]*x[2] - f[1]*x[0] + f[1]*x[2] + f[2]*x[0] - f[2]*x[1]))
    # X = sp.array([(x**2), x, sp.ones_like(x)]).T
    # a, b, c = l.solve(X,f.T)
#     sol = l.lstsq(X,f.T)
#     a,b,c = sol[0]
    # xmin = -b/(2*a)
    # This should never happen in the optimization case if the bracketing happened
    # correctly, but you've got to go to the minimum input point if the quadratic
    # is an "n" instead of a "u".
    # if sp.sign(a)==-1:
        # return min(x)
    return xmin

def GRatioOnce(x,f):
    return 1

def setBounds(ffile, xpoint, fpoint, stride,disp=0):
    if disp > 1.01: print('Bounding polynomial search')
    # Take 1 point on a function and the stride length, find three points
    # on the function that bracket the minimum.
    maxinner = 5
    maxouter = 50
    gamma = 0.5
    for bigi in range(maxouter):
        if not bigi == 0:
            xpoint = xq[-1]
            fpoint = fq[-1]
        bigi += 1
        sign = 1
        q = 0
        xq = []
        fq = []
        for d in range(maxinner):
            xq.append(xpoint + sign*(q)*stride)
            fq.append(ffile(xq[q]))
            if q==1 and fq[q]>fq[q-1]: # going uphill, reverse direction and rearrange
                sign = -1
                q = 2
                xq.append(xpoint + sign*(q)*stride)
                fq.append(ffile(xq[q]))
                temps = [xq[1], fq[1]]
                xq[1], fq[1] = [xq[0], fq[0]]
                xq[0], fq[0] = temps
            if q>1 and fq[q-2] >= fq[q-1] and fq[q] > fq[q-1]:
                xout = xq[-3:]
                fout = fq[-3:]
                if disp > 1.01: print('Bounding complete')
                return [xout, fout, 0]
            q += 1
        stride /= gamma
    return [None, None, 1]

def GRmin(ffile,x0,eps,istride,imax):
    # x0 is the starting location for x
    # eps is a tolerance for when to stop refining
    # stride is the initial stride length for bracketing the min\
    # tmul is a refining parameter (0<tmul<1) that determines how aggressively to refine
    # imax is the max number of refine iterations
    # dependencies:
    # setBounds() returns 3 points for x and f with min bracketed. Only 2 are used.
    # GROnce() takes three points and assoc. f evals and returns the min if it were a
    # quadratic
    iflag=0
    f0=ffile(x0)

    # generate init bracket
    xbounds, fbounds, maxiter = setBounds(ffile,x0,f0,istride)
    x = [xbounds[0], None, None, xbounds[-1]]
    f = [fbounds[0], None, None, fbounds[-1]]
    a = (3-sp.sqrt(5))/2
    tol=1
    for i in range(imax):
        if(i==0):
            x[1] = x[0]+a*(x[3] - x[0])
            x[2] = x[3]-a*(x[3] - x[0])
            f[1] = ffile(x[1])
            f[2] = ffile(x[2])
        if f[1] > f[2]: # min lies between 1 and 3
            # update left bound
            x[0], f[0] = x[1], f[1]
            x[1], f[1] = x[2], f[2]
            x[2] = x[3] - a*(x[3] - x[0])
            f[2] = ffile(x[2])
        else:
            x[3], f[3] = x[2], f[2]
            x[2], f[2] = x[1], f[1]
            x[1] = x[0] + a*(x[3] - x[0])
            f[1] = ffile(x[1])
        tol=sp.absolute(max(f)-min(f))/max([1e-8,sp.absolute(min(f))])
        minindex = sp.argmin(f)
        fstar = min(f)
        xstar = x[minindex]
        if(i==imax-1):
            print('warning max iteration reached in bracketer')
            iflag=1
        if tol < eps:
            break
    return [xstar, fstar, iflag]

def UPolyMin(func,g,x0,eps=1e-4,stride=0.1,tmul=0.5,imax=100,mode='GD',disp=0,pdisp=0,lineeps=1e-4,h=None,cmode='g'):
    # UPolyMin: Unconstrained Polynomial Line Search Minimization
    # Author: David Cunningham
    # # # # Inputs # # # # 
    # func  - function to be minimized. Must be of form: f(x):R^(nx)->R and accept only one (vector)
    #         argument
    # g     - gradient of function: Must be of form g(x):R^(nx)->R and accept only one (vector) arg
    # x0    - initial guess, x0 in R^(nx), array of floats
    # eps   - convergence criterion, float.
    #         When some error is smaller than this, you're done.
    # stride - Float. The step length the minimizer will travel to do the polynomial search
    #          on the first try.
    # tmul   - Float, 0<tmul<1. The multiplicative reduction factor for stride in case the poly-
    #          nomial search fails to bracket a minimum.
    # imax   - Int, the max iterations of the outer minimizer. Also currently controls the 
    #          max iterations of the line search.
    # mode   - String, chooses method for selecting line search direction.
    #          Possible values:
    #          - 'GD' : Gradient Descent
    #          - 'PR' : Polak Ribiere
    #          - 'FR' : Fletcher Reeves
    #          - 'BFGS' : Broyden-Fletcher-Goldfarb-Shanno
    # disp   - int, may be 0, 1, 2. Chooses the verbosity of the output for debugging.
    # pdisp  - int, may be 0, 1, 2. Chooses the verbosity of the line search output for debugging.
    # lineeps - convergence criterion for the inner line search, float.
    # cmode  - String, may be 'GD' or 'j'. Chooses what the stopping criterion is
    # # # # Outputs # # # #
    # x - list of guesses. The final element should be the converged guess, or the last guess
    # f - list of costs for guesses.
    # iflag - integer, will be 1 if the maximum iterations were hit anywhere. 

    x = []
    f = []
    iflag = 0
    x.append(x0)
    f.append(func(x0))
    for i in range(imax):
        if disp > 1.01: print('Computing direction')
        sh = dirf(mode,x,g,h)
        if disp > 1.01: print('Direction found')
        def line(t):
            xl = sh*t + x[-1]
            return func(xl)
        if disp > 1.01: print('Starting linesearch')
        amin,fstar,iflag = polymin(line, 0, lineeps, stride, tmul, imax, pdisp)
        if disp > 1.01: print('linesearch complete')
        x.append(x[i]+amin*sh)
        if disp > 0:
            print('Iteration:',i,'cost: ', fstar, 'Gradient norm: ', l.norm(dirf.ghist[-1]))
        f.append(fstar)
        if i==imax-1:
            print('warning max iteration in minimizer')
            iflag = 1
        if l.norm(dirf.ghist[-1]) < eps and cmode=='g':
            break
        if abs(f[-1])<eps and cmode=='J':
            break
    return [x, f, iflag]


def dirf(mode,x,g,sprev=None,h=None):
    # Keep a record of linesearch directions
    # Also keep a record of gradients
    if not hasattr(dirf, "hist"):
        dirf.hist = []
    if not hasattr(dirf, 'ghist'):
        dirf.ghist = []
    if mode=='GD':
        thisg = g(x[-1])
        dirf.hist.append(-thisg)
        dirf.ghist.append(thisg)
    elif mode=='FR':
        if len(x)==1:
            thisg = g(x[0])
            dirf.hist.append(-thisg)
            dirf.ghist.append(thisg)
        else:
            thisg = g(x[-1])
            dirf.ghist.append(thisg)
            gim1 = sp.array(dirf.ghist[-2])
            gi   = sp.array(dirf.ghist[-1])
            B    = (gi.T@gi)/(gim1.T@gim1)
            dirf.hist.append(-gi + B*dirf.hist[-1])

    elif mode=='PR':
        if len(x)==1:
            thisg = g(x[0])
            dirf.hist.append(-thisg)
            dirf.ghist.append(thisg)
        else:
            thisg = g(x[-1])
            dirf.ghist.append(thisg)
            gim1 = sp.array(dirf.ghist[-2])
            gi   = sp.array(dirf.ghist[-1])
            B    = (gi.T@(gi-gim1))/(gim1.T@gim1)
            dirf.hist.append(-gi + B*dirf.hist[-1])
    elif mode=='BFGS':
        thisg = g(x[-1])
        dirf.ghist.append(thisg)
        if not hasattr(dirf,"Qhist"):
            dirf.Qhist = [sp.eye(len(x[-1]))]
        if len(x)==1:
            dirf.Qhist = [sp.eye(len(x[-1]))]
        else:
            Qim1 = dirf.Qhist[-1]
            gim1 = sp.array(dirf.ghist[-2])
            gi   = sp.array(dirf.ghist[-1])
            p = x[-1] - x[-2]
            y = gi - gim1

            A = Qim1@sp.outer(y,p)
            tau = y@Qim1@y
            sig = p@y
            dQ  = (sig + tau)/tau**2 * sp.outer(p,p) - 1/sig*(A+A.T)
            dirf.Qhist.append(Qim1+dQ)
        dirf.hist.append(-dirf.Qhist[-1]@sp.array(dirf.ghist[-1]))
    else:
        raise Exception('Mode not implemented')
    return dirf.hist[-1]/l.norm(dirf.hist[-1])

def secantSolve(fun, x0, x1, eps=1e-4):
    # use secant method to rootsolve for lambda
    xs = [x0, x1]
    vals = [fun(x0), fun(x1)]
    tol = 1e5
    count = 0
    while tol > eps and count < 100:
        xs.append(xs[-1] - vals[-1]*(xs[-1] - xs[-2])/(vals[-1] - vals[-2]))
        vals.append(fun(xs[-1]))
        tol = abs(vals[-1])
        count +=1
    return xs[-1]


def matTens(M,T): #YES
    return np.einsum('ip,jpk -> jik',M,T)

def quad(M,T): #YES
    return np.einsum('pj,qk,piq -> jik',M,M,T)

def vecTens(v,T): #YES
    return np.tensordot(v,np.transpose(T,(1,2,0)),axes=1)

def adProp(x,Jac,STM,Hess=None,STT=None,order=1):
    # Takes square STM/STT returns vectorized STM/STT dot
    fx = Jac(x)
    STMdot = fx@STM
    dotvec = sp.reshape(STMdot,-1)
    # Note this tensor multiplication was checked against classmates and provided Matlab
    # (in Matlab)
    if order == 2:
        fxx = Hess(x)
        STTdot = mattens(fx,STT) + quad(STM,fxx)
        dotvec = sp.r_[dotvec, sp.reshape(STTdot,-1)]
    return dotvec
