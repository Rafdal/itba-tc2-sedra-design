import sympy as sp
from scipy import signal

def sympy_to_tf(H_s, s):
    # convert SymPy expression to transfer function
    num, den = sp.fraction(H_s)
    num_coeffs = sp.Poly(num, s).all_coeffs()
    den_coeffs = sp.Poly(den, s).all_coeffs()
    num_coeffs = [float(c) for c in num_coeffs]
    den_coeffs = [float(c) for c in den_coeffs]
    H_tf = signal.TransferFunction(num_coeffs, den_coeffs)
    return H_tf

def sympy_to_phase_gd(H_s, s, w):
    # convert SymPy expression to phase and group delay
    Hw = H_s.subs(s, sp.I*w)

    real, imag = sp.re(Hw), sp.im(Hw)

    phase = sp.atan2(imag, real).simplify()
    group_delay = -sp.diff(phase, w).simplify()
    return phase, group_delay

def sympy_format_polynomial(H_s, s):
    H_s = H_s.simplify()
    num, den = sp.fraction(H_s)
    num = sp.Poly(num, s)
    den = sp.Poly(den, s)
    return num/den