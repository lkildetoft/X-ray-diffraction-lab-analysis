import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
from scipy import signal
rc('text', usetex=True)

"""
quantitative analysis of silver data
"""

ag = pd.read_excel("data/Ag_xrd.xlsx", header=None) #read excel data file
ag_xdata = np.array(ag[0]) #split columns into arrays
ag_ydata = np.array(ag[1])

peaksag, _ = signal.find_peaks(ag_ydata, height = 45, distance=250) #find peaks
peaksag_clean = np.delete(peaksag, [0, 1]) #delete false peaks
fwhm = signal.peak_widths(ag_ydata, peaksag_clean, rel_height=0.5)

np.savetxt("2theta_peaks_ag.txt", ag_xdata[peaksag_clean]) #save 2theta values as txt

def scherrer(): #define scherrer function for task 7b

    global t
    global wl
    k = 0.94
    wl = 1.54e-10

    fwhmlist = fwhm[0]
    beta = sum([wide*((np.pi)/180) for wide in fwhmlist])/len([wide*((np.pi)/180) for wide in fwhmlist])

    theta = sum([(twotheta/2)*((np.pi)/180) for twotheta in ag_xdata[peaksag_clean]])/len([(twotheta/2)*((np.pi)/180) for twotheta in ag_xdata[peaksag_clean]])

    t = (k*wl)/(beta*np.cos(theta))

    return t

def miller(): #miller function for indexing reflections, returns rounded values of h^2 + k^2 + l^2

    global indices_list
    lattice_const = 4.086e-10

    indices_root = [(2*lattice_const*np.sin((twotheta/2)*((np.pi)/180)))/wl for twotheta in ag_xdata[peaksag_clean]]
    indices_list = [round(indices**2) for indices in indices_root]

    return indices_list

scherrer()
miller()

file = open("t.txt", "w+") #save calculated scherrer value to t.txt
file.write(str(t))
file.close()

file = open("miller.txt", "w+") #save calculated milller values to t.txt
file.write(str(indices_list))
file.close()

"""
mixture
"""
mix = pd.read_excel("data/mixture_xrd.xlsx", header=None) #read excel data file
mix_xdata = np.array(mix[0]) #split columns into arrays
mix_ydata = np.array(mix[1])

"""
aluminium oxide
"""
al = pd.read_excel("data/Al2O3_xrd.xlsx", header=None) #read excel data file
al_xdata = np.array(al[0]) #split columns into arrays
al_ydata = np.array(al[1])

"""
plotting
"""
plt.figure(1)
plt.plot(ag_xdata, ag_ydata)
plt.scatter(ag_xdata[peaksag_clean], ag_ydata[peaksag_clean], marker="x", color="orange")
plt.title("X-ray diffraction spectrum of silver sample, with peaks")
plt.xlabel(r"$2\theta$ [degrees]")
plt.ylabel("Intensity [arb. units]")
plt.savefig("silver.png", bboxinches="tight")
plt.show()

plt.figure(2)
plt.plot(ag_xdata, ag_ydata, label="Ag sample")
plt.plot(al_xdata, al_ydata, label="Al2O3 sample")
plt.plot(mix_xdata, mix_ydata, label="Mixture")
plt.legend()
plt.title("Overlayed X-ray diffraction patterns of three samples")
plt.xlabel(r"$2\theta$ [degrees]")
plt.ylabel("Intensity [arb. units]")
plt.savefig("overlay.png", bboxinches="tight")
plt.show()
