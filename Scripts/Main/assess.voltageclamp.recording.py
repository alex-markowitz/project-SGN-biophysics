### DESCRIPTION ###################################################################################
# SGN Recording EDA
### PREAMBLE #########################################################################

import pandas as pd
import pyabf
import matplotlib.pyplot as plt

abf = pyabf.ABF('/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/USC/PythonCoding/CellA004.abf');

plt.figure(figsize=(8, 5))
for sweepNumber in abf.sweepList:
    abf.setSweep(sweepNumber)
    i1, i2 = 0, int(abf.sampleRate * 1)  # plot part of the sweep
    dataX = abf.sweepX[i1:i2] + .025 * sweepNumber
    dataY = abf.sweepY[i1:i2] + 15 * sweepNumber
    plt.plot(dataX, dataY, color='C0', alpha=.5)

plt.gca().axis('off')  # hide axes to enhance floating effect
plt.show()