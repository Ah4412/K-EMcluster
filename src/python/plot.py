def plot(x, y, px, py):
    import matplotlib.pyplot as plt
    plt.scatter(x,y, s=150, c="green", alpha=0.5)
    plt.scatter(px,py, s=5, c="blue", alpha=0.5)
    plt.show()
    return 1