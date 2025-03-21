
import chem_analysis as ca


def main():
    ca.global_config.preferred_plot = "matplotlib"
    fig = ca.plot.signal(fid)
    fig = ca.plot.add_peaks(fid, bounds, fig=fig, labels=peak_labels, label_mode=1)
    fig.axes[0].set_xlim(4, 30)
    fig.axes[0].set_ylim(-10_000, 160_000)
    fig.show()


if __name__ == "__main__":
    main()
