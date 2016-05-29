def plot_xy(x, y):
    import matplotlib.pyplot as plt
    print "Printing from Python in plotStdVectors()"
    plt.scatter(x,y, s=150, c="green", alpha=0.5)
    plt.show()
    return 1