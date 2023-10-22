import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sp
from scipy import signal

class PlotBode:
    def __init__(self, figsize=(10, 6)):
        self.fig, self.ax1 = plt.subplots(figsize=figsize)
        self.ax2 = self.ax1.twinx()

        self.min_f = None
        self.max_f = None
        self.y1_limits = None
        self.y2_limits = None

    def __set_f_lim__(self, f):
        self.min_f = np.min([self.min_f, np.min(f)]) if self.min_f else np.min(f)
        self.max_f = np.max([self.max_f, np.max(f)]) if self.max_f else np.max(f)

    def __set_minimun_y_limits__(self):
        if self.y1_limits:
            ymin, ymax = self.ax1.get_ylim()
            if ymin > self.y1_limits[0]:
                self.ax1.set_ylim(bottom=self.y1_limits[0])
            if ymax < self.y1_limits[1]:
                self.ax1.set_ylim(top=self.y1_limits[1])
        if self.y2_limits:
            ymin, ymax = self.ax2.get_ylim()
            if ymin > self.y2_limits[0]:
                self.ax2.set_ylim(bottom=self.y2_limits[0])
            if ymax < self.y2_limits[1]:
                self.ax2.set_ylim(top=self.y2_limits[1])

    def plotTransfer(self, H_tf, f=None, w=None, points=1000):
        if w is None and f is None:
            raise ValueError('Either f or w must be specified')
        if w is None:
            w = 2*np.pi*f            

        x = w if (f is None) else f

        w, mag, phase = signal.bode(H_tf, w=w)

        self.plotSemilog1(x, mag, label='|H|_{{dB}}', color='red')
        self.plotSemilog2(x, phase, label='Fase [\\deg]', color='blue')
        self.y1_limits = [-40, 10]
        self.y2_limits = [-180, 180]

    def drawRectangle1(self, f1, f2, y1, y2, **kwargs):
        # self.ax1.axvspan(f1, f2, y1, y2, **kwargs)
        self.ax1.fill_between([f1, f2], y1, y2, **kwargs)

    def drawRectangle2(self, f1, f2, y1, y2, **kwargs):
        # self.ax2.axvspan(f1, f2, y1, y2, **kwargs)
        self.ax2.fill_between([f1, f2], y1, y2, **kwargs)

    def plotLinear1(self, f, y, **kwargs):
        self.ax1.plot(f, y, **kwargs)
        self.__set_f_lim__(f)

    def plotLinear2(self, f, y, **kwargs):
        self.ax2.plot(f, y, **kwargs)
        self.__set_f_lim__(f)

    def plotSemilog1(self, f, y, **kwargs):
        self.ax1.semilogx(f, y, **kwargs)
        self.__set_f_lim__(f)

    def plotSemilog2(self, f, y, **kwargs):
        self.ax2.semilogx(f, y, **kwargs)
        self.__set_f_lim__(f)
        
    def plotLoglog1(self, f, y, **kwargs):
        self.ax1.loglog(f, y, **kwargs)
        self.__set_f_lim__(f)

    def plotLoglog2(self, f, y, **kwargs):
        self.ax2.loglog(f, y, **kwargs)
        self.__set_f_lim__(f)

    def show(self, loc='best', min2loc=10, maj2loc=30, y1limits=None, gridcolor='gray', gridalpha=None):
        self.ax1.tick_params(axis='y', labelcolor='black')
        self.ax1.grid(True, which="both", ls="-", axis="both", color=gridcolor, alpha=gridalpha)

        self.ax2.tick_params(axis='y', labelcolor='black')
        self.ax2.grid(True, which="major", ls="-", color=gridcolor, alpha=gridalpha)

        self.ax1.xaxis.set_minor_locator(plt.LogLocator(base=10, subs='all', numticks=800))
        self.ax1.xaxis.set_major_locator(plt.LogLocator(base=10, numticks=200))

        # set ticks for ax2 y axis
        self.ax2.yaxis.set_major_locator(plt.MultipleLocator(maj2loc))
        self.ax2.yaxis.set_minor_locator(plt.MultipleLocator(min2loc))
        # self.ax1.yaxis.set_major_locator(plt.MultipleLocator(20))

        self.ax1.set_ylabel('Ganancia $|Z|_{{dB}}$')
        if self.ax2.get_lines():
            self.ax2.set_ylabel('Fase $[\degree]$')

        # set legend
        lines, labels = self.ax1.get_legend_handles_labels()
        lines2, labels2 = self.ax2.get_legend_handles_labels()
        self.ax2.legend(lines + lines2, labels + labels2, loc=loc)

        self.ax1.set_xlabel('Frecuencia $[Hz]$')

        # set axis limits
        self.ax1.set_xlim(left=self.min_f, right=self.max_f)
        self.ax1.set_ylim(y1limits)

        self.__set_minimun_y_limits__()
        plt.tight_layout()
        plt.show()

def plotSensTable(sensTable, labels, title="", figsize=(5,1)):
    plt.figure(figsize=figsize)
    plt.bar(labels, sensTable)
    plt.title(title)
    plt.grid(True, axis='y', alpha=0.5)
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(0.25))
    # plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.05))
    plt.axhline(y=0, color='k', alpha=0.6, lw=0.5)
    plt.show()

class PoleZeroPlotter:
    def __init__(self, transfer_function=None, figsize=(3, 3)):
        # create a new figure and axis object
        self.fig, self.ax = plt.subplots(figsize=figsize)
        # clear the axis
        self.ax.cla()
        self.zeros = []
        self.poles = []
        self.polesKwargs = {}
        self.zerosKwargs = {}

        # add zeros and poles from transfer function
        if transfer_function:
            for z in transfer_function.zeros:
                self.zeros.append(z)
            for p in transfer_function.poles:
                self.poles.append(p)

    def drawText(self, text, position, **kwargs):
        # if kwargs is empty, use default values
        if not kwargs:
            kwargs = {'fontsize': 12, 'color': 'black'}
        self.ax.text(position[0], position[1], text, **kwargs)
            
    def drawCircle(self, center, radius, **kwargs):
        # if kwargs is empty, use default values
        if not kwargs:
            kwargs = {'fill': False, 'color': 'black', 'ls': '--', 'alpha': 0.5}
        circle = plt.Circle(center, radius=radius, **kwargs)
        self.ax.add_artist(circle)
    
    def addPoles(self, poles, **kwargs):
        if isinstance(poles, np.ndarray):
            poles = poles.tolist()
        self.poles += poles

    def addZeros(self, zeros):
        if isinstance(zeros, np.ndarray):
            zeros = zeros.tolist()
        self.zeros += zeros
    
    def show(self, xMajLoc=25000, xMinLoc=5000, yMajLoc=10000, yMinLoc=5000, loc='best', aspect=None):
        zerosReal = [z.real for z in self.zeros]
        zerosImag = [z.imag for z in self.zeros]
        polesReal = [p.real for p in self.poles]
        polesImag = [p.imag for p in self.poles]

        # find max value for axis limits
        xLims = [0,0]
        yLims = [0,0]
        for z in self.zeros:
            xLims[0] = min(xLims[0], z.real) - 0.01
            xLims[1] = max(xLims[1], z.real) + 0.01
            yLims[0] = min(yLims[0], z.imag) - 0.01
            yLims[1] = max(yLims[1], z.imag) + 0.01
        for p in self.poles:
            xLims[0] = min(xLims[0], p.real) - 0.01
            xLims[1] = max(xLims[1], p.real) + 0.01
            yLims[0] = min(yLims[0], p.imag) - 0.01
            yLims[1] = max(yLims[1], p.imag) + 0.01
        for i in range(2):
            xLims[i] = xLims[i]*1.07
            yLims[i] = yLims[i]*1.07
        # xLims[0] += -1
        # xLims[1] += 0.5

        plt.grid(True, which="both", axis="both")
        self.ax.axhline(y=0, color='k', alpha=0.5)
        self.ax.axvline(x=0, color='k', alpha=0.5)

        # Set tick spacing
        self.ax.xaxis.set_major_locator(plt.MultipleLocator(xMajLoc))
        self.ax.xaxis.set_minor_locator(plt.MultipleLocator(xMinLoc))
        self.ax.yaxis.set_major_locator(plt.MultipleLocator(yMajLoc))
        self.ax.yaxis.set_minor_locator(plt.MultipleLocator(yMinLoc))

        # set axis limits
        self.ax.set_xlim(xLims[0], xLims[1])
        self.ax.set_ylim(yLims[0], yLims[1])

        # Fix aspect ratio
        if aspect:
            self.ax.set_aspect(aspect, 'box')

        # plot the zeros and poles
        self.ax.scatter(zerosReal, zerosImag, marker='o', color='blue', label='zeros')
        self.ax.scatter(polesReal, polesImag, marker='x', color='red', label='poles')
        
        # add legend and labels
        if loc:
            self.ax.legend(loc=loc)
        plt.tight_layout()

        # show the plot
        plt.show()

def fixPhaseJumps(phase):
    # Find the jumps in the zin_p plot
    jumps = np.where(np.abs(np.diff(phase)) > 180)[0] + 1

    # Add or subtract multiples of 360 to each point in the phase plot
    for i in range(len(jumps)):
        if phase[jumps[i]] > phase[jumps[i]-1]:
            phase[jumps[i]:] -= 360
        else:
            phase[jumps[i]:] += 360
    return phase


def save_pdf(filename, fig=None):
    """Save to @filename with a custom set of file formats.
    
    By default, this function takes to most recent figure,
    but a @fig can also be passed to this function as an argument.
    """
    if fig is None:
        plt.savefig("%s.%s"%(filename, "pdf"))
    else:
        fig.savefig("%s.%s"%(filename, "pdf"))