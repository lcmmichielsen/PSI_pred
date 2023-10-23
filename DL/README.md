# Deep learning models

Our deep learning models are based on the [Saluki](https://doi.org/10.1186/s13059-022-02811-x) model. Therefore this repository is forked from the original [Basenji](https://github.com/calico/basenji) repository, which contains the code for the Saluki model. We adapted the 'rnann.py' and 'dataset.py' scripts in the 'basenji' folder such that they can work with exons as input. We also added 'exon_*.py' scripts in the 'bin' folder to generate datasets, train the models, test the models, and do ISM. In the tutorials folder, we explain how to use these.

### Installation
Installation is similar to installing Basenji:

These tools were developed with Python3 and various scientific computing dependencies, which you can see and install via requirements.txt for pip and environment.yml for [Anaconda](https://www.continuum.io/downloads). For each case, we kept TensorFlow separate to allow you to choose the install method that works best for you. The codebase is compatible with the latest TensorFlow 2, but should also work with 1.15.

The easiest way to install the package is as follows:
- Create a new conda environment
- Install TensorFlow using the instructions on their [website](https://www.tensorflow.org/install)
- Clone this repo
- Run `pip setup.py develop`

To verify the install, launch python and run

```
    import basenji
```

### Training the models
See the [Tutorials](./tutorial) folder for information about training your own models. Example input data can be downloaded from Zenodo.

### Analyzing the results
See the [Figure reproducibility](../Figures) folder for information on how to analyze your results.
