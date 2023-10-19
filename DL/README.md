# Deep learning models

Our deep learning models are based on the [Saluki](https://doi.org/10.1186/s13059-022-02811-x) model. Therefore this repository is forked from the original [Basenji](https://github.com/calico/basenji) repository. We adapted the 'rnann.py' and 'dataset.py' scripts in the 'basenji' folder such that they can work with exons as input. We also added 'exon_*.py' scripts to generate datasets, train the models, test the models, and do ISM. In the tutorials folder, we explain how to use these.

### Installation
Installation is similar to installing Basenji:

These tools were developed with Python3 and various scientific computing dependencies, which you can see and install via requirements.txt for pip and environment.yml for [Anaconda](https://www.continuum.io/downloads). For each case, we kept TensorFlow separate to allow you to choose the install method that works best for you. The codebase is compatible with the latest TensorFlow 2, but should also work with 1.15.

Run the following to install dependencies and Basenji with Anaconda.
```
    conda env create -f environment.yml
    conda install tensorflow (or tensorflow-gpu)
    python setup.py develop --no-deps
```

Alternatively, if you want to guarantee working versions of each dependency, you can install via a fully pre-specified environment.
```
    conda env create -f prespecified.yml
    conda install tensorflow (or tensorflow-gpu)
    python setup.py develop --no-deps
```

Or the following to install dependencies and Basenji with pip and setuptools.
```
    python setup.py develop
    pip install tensorflow (or tensorflow-gpu)
```

Then we recommend setting the following environmental variables.
```
  export BASENJIDIR=~/code/Basenji
  export PATH=$BASENJIDIR/bin:$PATH
  export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH
```

To verify the install, launch python and run
```
    import basenji
```
