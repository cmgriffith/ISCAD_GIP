import numpy as np
from functions import *
import matplotlib
matplotlib.use('TkAgg')  # or 'TkAgg'
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# load parameters from .json
file = 'ISCAD_parameters.json'
params = loadjson(file)
# derived_params(params,file) # add derived values to params
globals().update(params) # update params in current python file 

def analytical_plots(params, Ls_s, mutuals, Rs_ACs, Ls_m_arr):
    # Rs_AC vs f
    plt.subplot(2,2,1)
    plt.plot(Rs_ACs[:, 0], Rs_ACs[:, 1]*1e3)
    plt.xlim([Rs_ACs[0,0], Rs_ACs[-1,0]])
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    plt.grid(color = 'blue', linestyle = '--', linewidth = 0.5)
    plt.grid(True)
    plt.xlabel('frequency')
    plt.ylabel('Slot resistance, mOhm')
    plt.title('Slot resistance vs Fundamental Frequency')
    # plt.show()

    # # L_sm vs Qs !!!!!!!!!!!
    # plt.subplot(2,2,2)
    # sumslots = np.arange(1,round(Qs/2/p + 1))
    # figure1 = plt.figure()
    # axes1 = figure1.add_subplot(1, 1, 1)
    # axes1.plot(sumslots, Ls_s*mutuals[1:],'d')
    # plt.xlabel("Qs Slot Number")
    # plt.ylabel("Mutual Inductance")
    # plt.show()
    
    # L_sm vs Qs !!!!!!!!!!!
    plt.subplot(2,2,2)
    sumslots = np.arange(1,round(Qs/2/p + 1))
    # figure1 = plt.figure()
    plt.plot(sumslots, Ls_s*mutuals[1:],'d--')
    plt.grid(color = 'blue', linestyle = '--', linewidth = 0.5)
    plt.grid(True)
    plt.xlabel("Qs Slot Number")
    plt.ylabel("Mutual Inductance")
    # plt.show() 

    # Plot bar chart of Ls_m vs Qs
    plt.subplot(2,2,3)
    plt.bar(range(0,Qs), Ls_m_arr)
    plt.xlabel("Qs (slot count)")
    plt.ylabel("L_m [H] main inductance")
    plt.title("Main inductance vs Q_s")
    plt.axhline(
        y=Ls_s,
        color='red',
        linestyle='--',
        label="L_s self inductance"
    )
    plt.legend()

    plt.suptitle("Analytical Plots")
    plt.show()


def plot_B(B_ag): # plot airgap flux
    print("plotting...")
    N = len(B_ag); # 361?
    Elect_angle = np.arange(N,dtype=int)*360/(N-1); # 0 to 360
    # matplotlib.use('TkAgg')
    figure1 = plt.figure()
    axes1 = figure1.add_subplot(1, 1, 1)
    axes1.set_xlim([0, 360])
    axes1.set_xticks([0, 60, 120, 180, 240, 300, 360])
    axes1.plot(Elect_angle, B_ag)
    # plot(Elect_angle,B_ag,Elect_angle,real(Bg_k(1)*exp((0:N-1)*2i*pi/(N-1))));
    # axes1.plot(Elect_angle, np.real(Bg_k[0] * np.exp(np.arange(0, N) * 2j * np.pi / (N - 1))))
    plt.xlabel("Angle [deg]")
    plt.ylabel("Flux density [T]")
    plt.show()

    plt.savefig("B airgap.png", transparent=None)


def plot_DFT(Bg_k):
    # fourier plot
    figure2 = plt.figure()
    axes2 = figure2.add_subplot(1, 1, 1)
    axes2.set_xlim([0, 14])
    axes2.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13])
    plt.bar(range(1,len(Bg_k)+1), np.abs(Bg_k)) # 0 to 12, want 1 to 13
    plt.xlabel("Harmonic Order")
    plt.ylabel("Magnitude [T]")
    plt.show()

    plt.savefig("B spectrum.png", transparent=None)


def plot_Aph(Aph):
    figure3 = plt.figure()
    axes3 = figure3.add_subplot(1, 1, 1)
    axes3.set_xlim([0, len(Aph)])
    axes3.set_xticks(np.arange(0,len(Aph),len(Aph)/10))
    print("plot_Aph")
    print(np.arange(0,len(Aph)))
    print(Aph)
    plt.bar(np.arange(0,len(Aph)), Aph)
    plt.xlabel("Qs Slot Number")
    plt.ylabel("Current [A]")
    plt.show()
    
    plt.savefig("Qs currents.png", transparent=None)


def plot_mutuals(Qs, mutuals):
    Qs = np.arange(0,Qs-1)
    figure4 = plt.figure()
    axes4 = figure4.add_subplot(1, 1, 1)
    print("plot_mutuals()")
    print(Qs)
    print(mutuals)
    axes4.plot(Qs, mutuals)
    plt.xlabel("Qs Slot Number")
    plt.ylabel("Current [A]")
    plt.show()

