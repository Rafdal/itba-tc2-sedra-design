import numpy as np
import sympy as sp


class RCHiPass:
    def __init__(self, w0, k):
        self.w0 = w0
        # check if k>1
        if k <= 1:
            raise ValueError('k must be greater than 1')
        self.k = k

    def setBaseComponents(self, C=None, R=None, Ra=None, Rb=None):
        # Just R or C can be set and just Ra or Rb can be set
        if (C is None and R is None) or (C is not None and R is not None):
            raise ValueError('C or R must be specified, not both or none')
        if (Ra is None and Rb is None) or (Ra is not None and Rb is not None):
            raise ValueError('Ra or Rb must be specified, not both or none')

        if C is None:
            self.R = R
            self.C = 1/(R*self.w0)
        else:
            self.C = C
            self.R = 1/(C*self.w0)
        if Ra is None:
            self.Rb = Rb
            self.Ra = (self.k - 1)*self.Rb
        else:
            self.Ra = Ra
            self.Rb = self.Ra/(self.k - 1)

    def printComponentValues(self):
        spiceStr = ".param"
        spiceStr += f' C = {self.C:.3e}'
        spiceStr += f' R = {self.R:.3e}'
        spiceStr += f' Ra = {self.Ra:.3e}'
        spiceStr += f' Rb2 = {self.Rb:.3e}'
        print(spiceStr)

    def printEngineeringFormatParams(self):
        # Engineering format uses M, k, m, u, n, p, f
        # M = 1e6, k = 1e3, m = 1e-3, u = 1e-6, n = 1e-9, p = 1e-12, f = 1e-15
        params = [self.C, self.R, self.Ra, self.Rb]
        labels = ['C', 'R', 'Ra', 'Rb2']

        for label, param in zip(labels, params):
            if param >= 1e6:
                print(f"{label} = {param/1e6:.2f}M")
            elif param >= 1e3:
                print(f"{label} = {param/1e3:.2f}k")
            elif param >= 1:
                print(f"{label} = {param:.2f} ")
            elif param >= 1e-3:
                print(f"{label} = {param*1e3:.2f}m")
            elif param >= 1e-6:
                print(f"{label} = {param*1e6:.2f}u")
            elif param >= 1e-9:
                print(f"{label} = {param*1e9:.2f}n")
            elif param >= 1e-12:
                print(f"{label} = {param*1e12:.2f}p")
            else:
                print(f"{label} = {param*1e15:.2f}f")

    def getTransferFunction(self, s=None):
        if s is None:
            s = sp.symbols('s')
        self.H_s = (self.k * s) /(s + self.w0)
        return self.H_s, s