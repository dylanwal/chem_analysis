# Chem Analysis

Your one-stop shop for analyzing chemistry data.

This package was developed for general chemistry data analysis, but special attention was paid to analyze large data sets 
(example: analyzing 1000 NMR at once).  

Design Philosophy:
* Handle large data sets 
  * Support analyzing 1000s of NMR at once. (like those generated from kinetic analysis)
* Modular 
  * Be able to turn on or off or switch out processing or analysis methods with minimal code change
* Explict
  * Alot of analytic software automatically perform data transforms and hide this from the user making it hard to 
  truly know what's going on with your data. Here everything needs to be called explicitly, but typical processing
  steps are suggested in several of the examples.


**Support data types**:
* IR
* NMR (Bruker, Spinsolve) - 1D only
* SEC (GPC)
* GC and LC/HPLC/UPLC
* Mass Spec. and GC-MS
* UV-Vis (coming soon)


## Installation
[pypi page](https://pypi.org/project/chem-analysis/)

`pip install chem_analysis`

### Dependencies
#### Required
* [numpy](https://github.com/numpy/numpy) (for math)
* [scipy](https://github.com/scipy/scipy) (for math)
* [tabulate](https://github.com/astanin/python-tabulate) (for pretty table outputs)
* [bigsmiles](https://github.com/dylanwal/BigSMILES) (for chemistry (molar weights, SMILES))
#### Optional
* [matplotlib](https://github.com/matplotlib/matplotlib) (for plotting)
* [plotly](https://github.com/plotly/plotly.py) (for plotting)
* [pyqtgraph](https://github.com/pyqtgraph/pyqtgraph) (for plotting)

## Capabilities
### Processing Methods:
* Baseline correction
* Peak Picking
* Resampling
* Translations
* Smoothing
* Phase correction (NMR)
* Referencing (NMR)
* Chromatogram Library Matching (GC, LC/HPLC/UPLC)
* Mass Spec Library Matching (MS, GC-MS)
* And more ...

### Analysis Methods:
* Peak picking
* Integration
* Peak fitting
* Deconvolution (peaks and MS spectra)
* Multi-component analysis (MCA)
* And more ...

## Plotting / GUI
* Matplotlib
  * Popular Python Plotting Library
* Plotly
  * Provides interactive plots (html)
  * Zoom in/out is a game changer
* PyQt
  * Uses GPU for rendering.
  * Good for lots of data!!

Plotting is hard and everyone has their preferences. The goal of plotting in this package if to 
provide quick and useful views of the data with useful defaults. 
It does not strive for perfect looks or full customization.  


## Quick Examples

See [Examples folder] for a large list.

```python
import pathlib
import plotly.graph_objs as go
import chem_analysis as ca


def main():
    cal_RI = ca.sec.ConventionalCalibration(lambda time: 10 ** (-0.6 * time + 10.644),
                                            mw_bounds=(160, 1_090_000), name="RI calibration")

    # loading data
    file_path = pathlib.Path(r"data//SEC.csv")
    sec_signal = ca.sec.SECSignal.from_csv(file_path)
    sec_signal.calibration = cal_RI

    # processing
    sec_signal.processor.add(
        ca.processing.baseline.ImprovedAsymmetricLeastSquared(lambda_=1e6, p=0.15)
    )

    # analysis
    peak = ca.analysis.peak_picking.find_peak_largest(sec_signal,
                                                      mask=ca.processing.weigths.Spans((10, 12.2), invert=True)
                                                      )
    result = chem_analysis.analysis.peaks.integration.rolling_ball(peak, n=45, min_height=0.05,
                                                                   n_points_with_pos_slope=1)

    # plotting
    fig = go.Figure(layout=ca.plotting.PlotlyConfig.plotly_layout())
    ca.plotting.signal(sec_signal, fig=fig)
    ca.plotting.calibration(cal_RI, fig=fig)
    ca.plotting.peaks(result, fig=fig)
    fig.data[0].line.color = "blue"  # customize colors
    fig.data[4].fillcolor = "gray"  # customize colors
    # fig.show()
    fig.write_html('figs/sec_data_analysis.png')

    # print results
    print(result.stats_table().to_str())


if __name__ == '__main__':
    main()


```

|   peak | area | ...** |   mw_d  |   mw_n  |
|--------|------|-------|--------|--------|
|      0 | 2.7  | ...   | 1.218  |   7465  |

**Only showing 4 of 28 stats calculated for SEC peak

![sec_data_analysis.png](https://github.com/dylanwal/chem_analysis/tree/develop/dev/sec_data_analysis.png)


## Important philosophy for working with the code base
* Both `signal` and `x, y` options are available for all processing methods.
  * `Signal` are merely a convince object that can be bringing important information into analysis and plotting
* 'Spans' are 'x' values and 'Slices' are index values. The default is to prefer 'Spans'

## Contributing

Let's be honest, that are bugs in the code! Pointing them out are much appreciated, and solutions are double appreciated! 

Got a new data type, new algorithm, analysis, etc.? New contributions are welcomed!  
Best practice is to open an issue with your idea, and I will let you know if it
is a good fit for the project. If you are interested in helping code the addition please mention that as well. 

