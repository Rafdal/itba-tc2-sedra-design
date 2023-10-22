import sympy as sp
import numpy as np

# Clase de Python para diseñar una celda Sedra con Ceros de Transmisión y Polos con w0 y Q
# Diseñada para polos con Q < 5, en caso contrario se deberia utilizar otro enfoque

class SedraHPB:
    def __init__(self, w0, Q, wz):
        self.w0_val = w0
        self.Q_val = Q
        self.wz_val = wz

    def setDesignParams(self, Q0, n2):
        self.n2_lim = n2
        # Check if Q0 is valid (Q0 < Q) with a margin of 2%
        # if Q0*1.02 >= self.Q_val:
        #     raise ValueError('Q0 must be less than Q')
        self.Q0_val = Q0

    def setBaseComponentValues(self, Rb, C):
        self.Rb = Rb
        self.C = C

    def computeDesignParams(self):
        Q0 = self.Q0_val
        Q = self.Q_val
        n2 = self.n2_lim
        wz = self.wz_val
        w0 = self.w0_val
        self.K = (1/(2*Q0**2))*(1 - Q0/Q) + 1
        self.k1 = (n2*(wz/w0)**2)/(1 - Q0/Q)
        self.n = self.k1*(1 - Q0/(self.K*Q))
        self.m = self.k1*((self.K-1)/self.K)*(1 + (2*Q0**2)*(w0/wz)**2)

    def computeComponentValues(self):
        k = self.k1
        K = self.K
        n = self.n
        m = self.m
        Q0 = self.Q0_val
        Rb = self.Rb
        C = self.C
        Gb = (1/Rb)
        G = (self.w0_val*self.C)/(2*Q0)
        
        self.Gb = Gb
        self.Ga1 = (1-k)*(K-1)*Gb
        self.Ga2 = k*(K-1)*Gb
        self.G1 = 4*G*Q0**2
        self.G42 = n*G
        self.G41 = (1-n)*G
        self.C3 = C
        self.C22 = m*C
        self.C21 = (1-m)*C

        # Resistor values
        self.Ra1 = 1/self.Ga1
        self.Ra2 = 1/self.Ga2
        self.R1 = 1/self.G1
        self.R42 = 1/self.G42
        self.R41 = 1/self.G41
        # Rb is already defined        

    def printEngineeringFormatParams(self):
        # Engineering format uses M, k, m, u, n, p, f
        # M = 1e6, k = 1e3, m = 1e-3, u = 1e-6, n = 1e-9, p = 1e-12, f = 1e-15
        params = [self.C3, self.C22, self.C21, self.Ra1, self.Ra2, self.R1, self.R42, self.R41, self.Rb]
        labels = ['C3', 'C22', 'C21', 'Ra1', 'Ra2', 'R1', 'R42', 'R41', 'Rb']

        for label, param in zip(labels, params):
            if param >= 1e6:
                print(f"{label} = {param/1e6:.3f} M")
            elif param >= 1e3:
                print(f"{label} = {param/1e3:.3f} K")
            elif param >= 1:
                print(f"{label} = {param:.3f} ")
            elif param >= 1e-3:
                print(f"{label} = {param*1e3:.3f} m")
            elif param >= 1e-6:
                print(f"{label} = {param*1e6:.3f} u")
            elif param >= 1e-9:
                print(f"{label} = {param*1e9:.3f} n")
            elif param >= 1e-12:
                print(f"{label} = {param*1e12:.3f} p")
            else:
                print(f"{label} = {param*1e15:.3f} f")

    def printSpiceParams(self):
        spiceStr = ".param"
        sortList = [
            [f" C3={self.C3:.3e}",  self.C3],
            [f" C22={self.C22:.3e}",  self.C22],
            [f" C21={self.C21:.3e}",  self.C21],
            [f" Ra1={self.Ra1:.3e}",  self.Ra1],
            [f" Ra2={self.Ra2:.3e}",  self.Ra2],
            [f" R1={self.R1:.3e}",  self.R1],
            [f" R42={self.R42:.3e}",  self.R42],
            [f" R41={self.R41:.3e}",  self.R41],
            [f" Rb={self.Rb:.3e}",  self.Rb],
        ]
        sortList.sort(key=lambda x: x[1])

        for el in sortList:
            spiceStr += el[0]
        print(spiceStr)

    def packComponentValues(self):
        self.computeDesignParams()
        self.computeComponentValues()
        self.values = {
            'Gb': self.Gb,
            'Ga1': self.Ga1,
            'Ga2': self.Ga2,
            'G1': self.G1,
            'G42': self.G42,
            'G41': self.G41,
            'C3': self.C3,
            'C22': self.C22,
            'C21': self.C21,
            'Ra1': self.Ra1,
            'Ra2': self.Ra2,
            'R1': self.R1,
            'R42': self.R42,
            'R41': self.R41,
            'Rb': self.Rb
        }
        return self.values
    
    def getTransferFunction(self, s=None):
        if s is None:
            s = sp.symbols('s')
        Ga1 = self.Ga1
        Ga2 = self.Ga2
        G1 = self.G1
        G42 = self.G42
        G41 = self.G41
        Gb = self.Gb
        C3 = self.C3
        C22 = self.C22
        C21 = self.C21
        
        # Denominador (Polos)
        self.w02 = G1*(G41 + G42)/(C3*(C21 + C22))
        self.w0_Q = (G41 + G42)*(1/(C21 + C22) + 1/C3) - (G1/(C21 + C22))*(Ga1 + Ga2)/Gb
        self.n2 = ((Ga1 + Ga2 + Gb)/Gb)*((C22)/(C21 + C22)) - Ga2/Gb
        self.n1 = (Ga1/Gb + Ga2/Gb + 1)*G42*(1/(C21 + C22) + 1/C3) - (Ga2/Gb)*( G1/(C21 + C22) + (G41 + G42)*(1/(C21 + C22) + 1/C3) )
        self.n0 = (G1*(G41 + G42)/(C3*(C21 + C22)))*((G42/(G41 + G42))*(Ga1/Gb + Ga2/Gb + 1) - Ga2/Gb )

        self.num_coeffs = [self.n2, self.n1, self.n0]
        self.den_coeffs = [1, self.w0_Q, self.w02]

        self.num = sp.Poly(self.num_coeffs, s)
        self.den = sp.Poly(self.den_coeffs, s)
        self.H = self.num/self.den
        return self.H, s
    
    def calculateSensTables(self, ponderar=True):
        k = self.k1
        K = self.K
        n = self.n
        m = self.m
        Q_0 = self.Q0_val
        Q = self.Q_val
        wz = self.wz_val
        w0 = self.w0_val
        n2 = self.n2_lim
        Ra1 = self.Ra1
        Ra2 = self.Ra2
        R1 = self.R1
        R42 = self.R42
        R41 = self.R41
        Rb = self.Rb
        C3 = self.C3
        C22 = self.C22
        C21 = self.C21

        self.w0_table = [-Q*n2*(2*Q_0**2*w0**2 + wz**2)/(2*w0**2*(2*Q*Q_0**2 + Q - Q_0)),
        Q*wz**2*n2*(-2*Q_0**2 - 1)/(2*w0**2*(2*Q*Q_0**2 + Q - Q_0)),
        -1/2,
        (-2*Q*Q_0**2*w0**2 + 2*Q*Q_0**2*wz**2*n2 - Q*w0**2 + Q*wz**2*n2 + Q_0*w0**2)/(2*w0**2*(2*Q*Q_0**2 + Q - Q_0)),
        -1/2,
        (Q*n2*(2*Q_0**2*w0**2 + wz**2) - w0**2*(2*Q*Q_0**2 + Q - Q_0))/(2*w0**2*(2*Q*Q_0**2 + Q - Q_0)),]
        # self.w0_table = [-1/2, -1/2, -m/2, -n/2, n/2 - 1/2, m/2 - 1/2]
        self.Q_table = [2*Q_0**2*k*(K - 1)/(2*K*Q_0**2 - 2*Q_0**2 - 1), Q_0**2*m*(K - 1)/(2*K*Q_0**2 - 2*Q_0**2 - 1), n*(-2*K*Q_0**2 + 2*Q_0**2 - 1)/(2*(2*K*Q_0**2 - 2*Q_0**2 - 1)), (K*Q_0**2*n - K*Q_0**2 - Q_0**2*n + Q_0**2 + n/2 - 1/2)/(2*K*Q_0**2 - 2*Q_0**2 - 1), (K*Q_0**2 - Q_0**2 + 1/2)/(2*K*Q_0**2 - 2*Q_0**2 - 1), 2*Q_0**2*(1 - K)/(2*K*Q_0**2 - 2*Q_0**2 - 1), 2*Q_0**2*(-K*k + K + k - 1)/(2*K*Q_0**2 - 2*Q_0**2 - 1), Q_0**2*(1 - K)/(2*K*Q_0**2 - 2*Q_0**2 - 1), Q_0**2*(-K*m + K + m - 1)/(2*K*Q_0**2 - 2*Q_0**2 - 1)]
        self.wz_table = [(m*(k*(K - 1)*(n - 1) - n*(K - 1)*(k - 1) + n) + n*(-k*(K - 1)*(m - 1) + m*(K - 1)*(k - 1) - m))*((K - 1)*(k - 1) - 1)/(2*(k*(K - 1)*(m - 1) - m*(K - 1)*(k - 1) + m)*(k*(K - 1)*(n - 1) - n*(K - 1)*(k - 1) + n)), n*(-K*k + K + k)/(2*(K*k - K*n - k)), m*(-K*k + K + k)/(2*(K*k - K*m - k)), k*(K*n - K - n + 1)/(2*(K*k - K*n - k)), -1/2, k*(-K*m + K*n + m - n)/(2*(K**2*k**2 - K**2*k*m - K**2*k*n + K**2*m*n - 2*K*k**2 + K*k*m + K*k*n + k**2)), (-(k*(K - 1)*(m - 1) + m)*(k*(K - 1)*(n - 1) - n*(K - 1)*(k - 1) + n) + (k*(K - 1)*(n - 1) + n)*(k*(K - 1)*(m - 1) - m*(K - 1)*(k - 1) + m))/(2*(k*(K - 1)*(m - 1) - m*(K - 1)*(k - 1) + m)*(k*(K - 1)*(n - 1) - n*(K - 1)*(k - 1) + n)), -1/2, k*(K*m - K - m + 1)/(2*(K*k - K*m - k))]

        self.w0_labels = ['C3', 'R1', 'C22', 'R42', 'R41', 'C21']
        self.Q_labels = ['Ra2', 'C22', 'R42', 'R41', 'R1', 'Rb', 'Ra1', 'C3', 'C21']
        self.wz_labels = ['Ra2', 'R42', 'C22', 'R41', 'R1', 'Rb', 'Ra1', 'C3', 'C21']

        # Calculate square sum of sensitivities
        self.w0_sen_sum = 0
        self.Q_sen_sum = 0
        self.wz_sen_sum = 0

        Rtol = 0.01
        Ctol = 0.1

        for i, sen in enumerate(self.w0_table):
            if ponderar:
                if 'R' in self.w0_labels[i]:
                    sen = sen*Rtol
                elif 'C' in self.w0_labels[i]:
                    sen = sen*Ctol
                sen2 = np.abs(sen)
                self.w0_table[i] = sen2
                self.w0_sen_sum += sen2
            else:
                self.w0_sen_sum += np.abs(sen)
        for i, sen in enumerate(self.Q_table):
            if ponderar:
                if 'R' in self.Q_labels[i]:
                    sen = sen*Rtol
                elif 'C' in self.Q_labels[i]:
                    sen = sen*Ctol
                sen2 = np.abs(sen)
                self.Q_table[i] = sen2
                self.Q_sen_sum += sen2
            else:
                self.Q_sen_sum += np.abs(sen)
        for i, sen in enumerate(self.wz_table):
            if ponderar:
                if 'R' in self.wz_labels[i]:
                    sen = sen*Rtol
                elif 'C' in self.wz_labels[i]:
                    sen = sen*Ctol
                sen2 = np.abs(sen)
                self.wz_table[i] = sen2
                self.wz_sen_sum += sen2
            else:
                self.wz_sen_sum += np.abs(sen)