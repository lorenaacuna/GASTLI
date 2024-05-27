

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt



class water_curves:
    def __init__(self, path_to_file):
        """
        Object that contains plots water phase curves for phase diagram


        Args:
            path_to_file: path to input data directory

        """


        # Arguments of __init__
        self.path_to_file = path_to_file

        data = pd.read_csv(path_to_file+'Input/Plotting/water_curves.dat', sep='\s+')

        self.P1 = data['P(1)']
        self.P2 = data['P(2)']
        self.P3 = data['P(3)']
        self.P4 = data['P(4)']
        self.P5 = data['P(5)']
        self.P6 = data['P(6)']
        self.P7 = data['P(7)']
        self.T = data['T']

        self.points = pd.read_csv(path_to_file+'Input/Plotting/water_tpts.dat', sep='\s+')

        #############################
        #############################
        # END Reading in files
        #############################
        #############################

    def plot_water_curves(self,ax):
        r"""Function that plots the lines of water phase diagram curves

        Args:
            ax:
                matplotlib axes class object where diagram is added

        Returns:
            ax:
                matplotlib axes class with final water diagram

        """

        ax.plot(self.T[0:27315], self.P1[0:27315], '-k')
        ax.plot(self.T[25116:27316], self.P2[25116:27316], '-k')
        ax.plot(self.T[25116:25615], self.P3[25116:25615], '-k')
        ax.plot(self.T[25616:27330], self.P4[25616:27330], '-k')
        ax.plot(self.T[27330:35499], self.P5[27330:35499], '-k')
        ax.plot(self.T[35499:-1], self.P6[35499:-1], '-k')
        ax.plot(self.T[27315:64708], self.P7[27315:64708], '-k')
        ax.plot([self.T[len(self.T) - 1], 1700], [self.P6[len(self.P6) - 1], 1e12], '--k')

        ax.plot([0.0, 2.511650e2], [2.085660e8, 2.085660e8], '-', color='tab:grey')
        ax.plot([3.550000E+02, 1.000000E+02], [2.216000E+09, 1.000000E+12], '-', color='tab:grey')
        ax.plot([6.470960E+02, 6.470960E+02], [2.206400E+07, 1.167296E+10], '-', color='tab:grey')
        ax.plot([6.470960E+02, 20000], [2.206400E+07, 2.206400E+07], '-', color='tab:grey')

        ax.plot(self.points['Ttp'], self.points['Ptp'], '.k')

        self.ax = ax

