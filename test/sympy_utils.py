import sympy as sym
from math import prod
def xi(i,name='x'):
    return sym.symbols(''.join([name,'(',i.__repr__(),')']))

def packed_lin_ind(l=7,si=0,lis=0,oind=1,name='x'):
    # Build list of tuples-> linear index conversions
    # ptups takes a (i,j) or (i,j,k) tuple and converts it to 
    # the linear index of a packed state with 0 and 1 appended to the end
    lr = l+si
    fullrange = range(si,lr)
    skipcolumn = range(si,lr-1)
    skiprow = range(si,lr-1)
    pind = lis
    ptups = {}
    inarray = [xi(i+1,name=name) for i in range(l)]
    for i in fullrange:
        pind = pind+1
        oind = oind+1
    for j in skipcolumn:
        for i in range(si,lr):
            ptups[(i,j)] = pind
            inarray.append(xi(oind,name=name))
            oind = oind+1
            pind = pind+1
    for k in skipcolumn: # Skip the last column of each page in the stt
        for j in range(k,lr-1): # lower triangular, skipping the last row
            for i in fullrange:
                   # Assign same index to both sides of a diagonal
                   ptups[(i,j,k)] = pind
                   ptups[(i,k,j)] = pind
                   inarray.append(xi(oind,name=name))
                   pind = pind+1
                   oind = oind+1

    zind = pind
    inarray.append(0)
    uind = pind+1
    inarray.append(1)

    for i in range(si,lr):
        ptups[(i,l-1+si)] = zind

    ptups[(l-1+si,l-1+si)] = uind

    # Outer column and row of STT are zero
    for i in range(si,lr):
        for j in range(si,lr):
            ptups[(i,j,l-1+si)] = zind
            ptups[(i,l-1+si,j)] = zind
    return (inarray,ptups)

def mattenstolindict(l=7,order=2,si=0,lis=0):
    matlinind = lis
    tenslinind = lis
    lr = l+si
    fullrange = range(si,lr)
    cm2 = {}
    cm3 = {}
    # cm2 accepts a (j,k) tuple as key and returns the index
    # of the same element in a column-major linear array
    # cm3 accepts a (i,j,k) tuple as key and returns the index
    # of the same element in a column-major linear array
    for k in fullrange:
        for j in fullrange:
            # Make a list of matrix indices in cmaj order
            cm2[(j+si,k+si)] = matlinind
            matlinind = matlinind + 1
            for i in fullrange:
                # list of tensor indices in cmaj order
                cm3[(i+si,j+si,k+si)] = tenslinind
                tenslinind = tenslinind + 1
    return (cm2, cm3)

def makeeoms(xprime,statevector,Amatlist,Hesslist,si=0,lis=0):
    # Outer iterator (alphabetical order reflects array index position)
    # Loops are ordered for column major linear array storage
    # For ease of use in Fortran
    matlinind = lis
    tenslinind = lis
    inarray, ptups  = packed_lin_ind(len(statevector),si=si,lis=lis)
    cm2,cm3 = mattenstolindict(len(statevector),order=2,si=si,lis=lis)
    eoms = [*xprime]

    # make sure we are using lower triangular stt/hessian elements only
    for k,_ in enumerate(statevector):
        for j,_ in enumerate(statevector):
            for i,_ in enumerate(statevector):
                if k>j: #upper triangular
                    Hesslist[cm3[(i+si,j+si,k+si)]] = \
                    Hesslist[cm3[(i+si,k+si,j+si)]]


    # Construct the STM EOMs
    for a in range(len(statevector)-1):
        for i in range(len(statevector)):
            terma = sum([Amatlist[cm2[(i+si,al+si)]]*\
                         inarray[ptups[(al+si,a+si)]]
                         for al in range(len(statevector))])
            eoms.append(terma)

    # Construct the STT EOMs
    for b in range(len(statevector)-1):
        for a in range(b,len(statevector)-1):
            for i in range(len(statevector)):
                terma = sum([Amatlist[cm2[(i+si,al+si)]]*\
                             inarray[ptups[(al+si,a+si,b+si)]]
                             for al in range(len(statevector))])
                termb = sum([sum([Hesslist[cm3[(i+si,al+si,be+si)]]*\
                                  inarray[ptups[(al+si,a+si)]]*\
                                  inarray[ptups[(be+si,b+si)]]
                             for al in range(len(statevector))])
                             for be in range(len(statevector))])
                eoms.append(terma+termb)

    return inarray,eoms

def makejh(xprime,statevector,H=None,order=2,pack=False,symm=True,simple=False):
    if order >2:
        print("order <=2 ... it says hess right there in the name")
        return
    if order >=1:
        Amatrix = []
        if H is not None:
            Hamgrad = []
    if order ==2:
        hess = []
        if H is not None:
            Hamhess = []
    # Outer iterator (alphabetical order reflects array index position)
    # Loops are ordered for column major linear array storage
    # For ease of use in Fortran
    for k, dxk in enumerate(statevector):
        # Hamiltonian part
        if H is not None:
            Hamgrad.append(H.diff(dxk))
        for j, dxj in enumerate(statevector):
            # Hamiltonian part
            if H is not None and order==2:
                Hamhess.append(H.diff(dxj,dxk))
            # Don't do anything if packed array
            # Note this continues over hessian block too
            if (pack \
                and j==len(statevector)):
                    continue
            # Amatrix part
            Amatrix.append(xprime[j].diff(dxk))

            # EOM part
            if order == 2:

                # Don't add to hessian if packing conditions are met
                for i, dxi in enumerate(statevector):
                    if (pack \
                        and k>j \
                        and k<len(statevector)):
                        continue
                    # If symmetry is on, make sure you use the symmetric element
                    elif symm and k<j:
                        # Lower triangular
                        hess.append(xprime[i].diff(dxk,dxj))
                    else:
                        # Upper triangular
                        hess.append(xprime[i].diff(dxj,dxk))

    returns = [Amatrix]
    if order == 2: returns.append(hess)
    if H is not None:
        returns.append(Hamgrad)
        if order == 2: returns.append(Hamhess)
    if simple:
        for r in returns:
            for obj in r:
                for o in obj:
                    o = o.simplify()
    return tuple(returns)

def make_fort(L,name="Default",shape=None,odesub=False,invars="x",prec='8'):
    from sympy.printing import fcode
    from sympy.simplify.cse_main import opt_cse
    print("")
    print("! BEGIN AUTOCODE OUTPUT FOR", name.upper())
    if shape is None:
        shape = [len(L)]
    optexpr = opt_cse(L,order="none")
    cses, mat = sym.cse(sym.Matrix(L).subs(optexpr.items()))
    mat = mat[0]
    dimstring = ','.join([str(n) for n in shape])
    csevarsep = 8
    csevarstring1 = "    real(" + prec + ")             :: &"

    restofcsevars = []
    for i in range(0,len(cses),csevarsep):
        restofcsevars.append('                          & '\
                             +', '.join([str(exp[0]) for exp in cses[i:i+csevarsep]]))

    if odesub:
        headertext = ''.join(["subroutine ",name,"(t,x,res,iflg)"])
        footertext = ''.join(["end subroutine ",name])
        result = "    real(" + prec + "), dimension("+ dimstring + "), intent(out) :: res"
    else:
        headertext = ''.join(["function ",name,"(",','.join(invars),") result(res)"])
        footertext = ''.join(["end function ",name])
        result = "    real(" + prec + "), dimension("+ dimstring + ") :: res"

    if len(shape)!=1:
        interstring = "    real(" + prec + "), dimension(" + str(prod(list(shape))) +") :: inter"
        assignvar = "inter"
    else:
        interstring = ""
        assignvar = "res"

    print(headertext)
    print("    implicit none")
    if odesub:
        print("    real(" + prec + "), intent(in) :: x(:),t")
        print("    integer,  intent(in) :: iflg")
    else:
        print("    real(" + prec + "), intent(in) :: " + '(:), '.join(invars),"(:)")
    print(csevarstring1)
    print(*restofcsevars,sep=', & \n')
    print(interstring)
    print(result)
    print("")
    for ce in cses:
        print("   ", \
              str(ce[0]),\
              " = ", \
              str(fcode(ce[1],standard=2003,source_format='free')).\
              replace("d0","_"+prec))
    print("")
    for i,m in enumerate(mat):
        print(''.join(['    ',assignvar,'(',str(i+1),') = ']),m)

    if assignvar=="inter":
        print(''.join(["    res = reshape(inter,[",*dimstring,"])"]))

    print(footertext)
    print("! END AUTOCODE OUTPUT FOR", name.upper())
