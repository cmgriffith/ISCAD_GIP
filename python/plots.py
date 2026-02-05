import numpy as np
import matplotlib.pyplot as plt

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
    plt.bar(np.arange(0,len(Aph)), Aph)
    plt.xlabel("Qs Slot Number")
    plt.ylabel("Current [A]")
    plt.show()
    
    plt.savefig("Qs currents.png", transparent=None)


def plot_mutuals(Qs, mutuals):
    Qs = np.arange(0,Qs-1)
    figure4 = plt.figure()
    axes4 = figure4.add_subplot(1, 1, 1)
    axes4.plot(Qs, mutuals)
    plt.xlabel("Qs Slot Number")
    plt.ylabel("Current [A]")
    plt.show()