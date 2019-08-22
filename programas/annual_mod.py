import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import math as math
import pandas as pd

# Labels with Latex:
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def main(save_bool = True):
    # Compute the velocity modulation:
    t = np.asarray([i for i in range(0, 366)])
    vE = 232 + 15*np.cos(2*np.pi* (t - 152.5)/365.25)

    # Define the dates:
    startday = pd.datetime(2019, 1, 1)
    days = pd.date_range(startday, periods=366, freq='D')

    # Plot:
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    xaxis = ax.get_xaxis()
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    plt.plot(days, vE)
    plt.xlabel('Months', fontsize=16)
    plt.ylabel('$V_E\ (m\cdot s^{-1})$', fontsize=16)
    plt.xticks(rotation=45, fontsize=12)
    plt.tight_layout()
    
    # Saving:
    if save_bool:
        plt.savefig('/afs/ciemat.es/user/t/tomas/CIEMAT/Dark_Matter/programas/figures/Annual_Mod.png')

    plt.show()
    
main()