import matplotlib.pyplot as plt


def display():
    hist = plt.hist2d(X, Y, bins=[23, 36])
    plt.colorbar(hist[3])
