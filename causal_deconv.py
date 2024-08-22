from numpy import poly1d
from sympy import Poly, symbols, simplify, O, fraction
from scipy.signal import welch, convolve
import statsmodels.api as sm

# thanks, chatGPT, for writing this helper function 
# that got the reverse polynomial ordering correct!
def numpy_to_sympy(numpy_poly):
    # Get the coefficients and degree of the NumPy polynomial
    coefficients = numpy_poly[::-1]
    degree = len(coefficients) - 1

    # Create a list of SymPy symbols
    x = symbols('x')

    # Create the SymPy polynomial using the coefficients and symbols
    sympy_poly = Poly(sum(coefficients[i] * x**i for i in range(degree + 1)), x, domain = 'CC')

    return sympy_poly

# This is where most of the heavy lifting is: 
# compute the ARIMA spectrum and Q(z) and Q(1/z) polynomials
def get_Qzs(fluor, p=1, d=1, q=1):

    # mainly working with z, with x as a dummy variable
    x, z = symbols("x z")
    
    # ARIMA(p,d,q) model
    # this means we're working with differenced data, so will need to integrate later
    model = sm.tsa.statespace.SARIMAX(fluor, order=(p, d, q)).fit()
    
    # ARMA coefficients correspond to polynomials in 1/z (forward time shifts)
    # have to turn the model coeff lists into polynomials in sympy
    # 
    b = numpy_to_sympy(model.polynomial_ar[::-1]).subs(x, 1/z)
    a = numpy_to_sympy(model.polynomial_ma[::-1]).subs(x, 1/z)
    Qz = b/a
    
    # maybe overkill on terms
    Qz_series= Qz.series(n =10).removeO()
    binv = numpy_to_sympy(model.polynomial_ar[::-1]).subs(x, z)
    ainv = numpy_to_sympy(model.polynomial_ma[::-1]).subs(x, z)

    Qzinv_series = (binv/ainv).series(n=10).removeO()

    variance = float(model.summary().tables[1].data[-1][1])
    return(b,a, binv, ainv, Qz_series, Qzinv_series, variance)

def expand(expr):
    return expr.series(n=10).removeO()

# remove all terms with exponent > 0
def causal_terms(expr):
    z = symbols("z")
    return (expr+O(z)).removeO()


# now build the filter
def build_cwf(fluor, T, falltime, W):
    
    # here is z, as a sympy symbol
    z = symbols("z")
    
    # some actual work is done here
    b,a, binv, ainv, Qz_series, Qzinv_series, variance = get_Qzs(fluor, p=2, d=1, q=2)
    
    # inverse of exponential in z domain (scaled based on sampling period)
    Ginv = (1-np.exp(-T/falltime)*1/z)

    # build filter using the spectral factorization
    filt_unsimplified = a/(b*variance)*causal_terms(expand(Ginv*variance*b/a)-expand(Ginv*W*ainv/binv))
    
    # clean up the algebra, make it a nice fraction in 1/z
    num, den = fraction(simplify(filt_unsimplified))
    z_filter = num/den
    
    # for my next trick i will turn a rational function into a filter
    # i intone the magic words "z is z, in whichsoever library it may be"
    # highly unsafe in general; do not try this at home
    
    from audiolazy import z
    H = eval(str(z_filter))
    
    return H

# causal Wiener deconvolution
def deconvolve_causal(trace, time, falltime, W = None):
    import matplotlib.pyplot as plt
    exp_fall = np.exp(-(time-time[0])/falltime)
    f, psd = welch(trace)
    if W is None:
        W = np.min(psd)
    T = time[1]-time[0]
    print(f'Deconvolving (causal method; |Wn| = {W})...')
    H = build_cwf(trace, T, falltime, W)
    deconvolved = stable_quadr(np.array(list(H(np.diff(trace)))))
    Hinv = 1/H
    reconst = stable_quadr(np.array(list(Hinv(np.diff(deconvolved, prepend = 0)))))

    return deconvolved, reconst, H, W

# mean-zero reintegration
def stable_quadr(x):
    return np.cumsum(x-np.mean(x))
