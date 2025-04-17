# Overview

This directory contains code and artifacts for training a neural network to
approximate droplet characteristics.  This encompasses data generation, training
via PyTorch, and creation of a Fortran module allows inference with model to
use it in simulations.

Contents of this directory include:

* `bin/`: Collection of scripts to run from the command line
* `data/`: Where training data are generated
* `models/`: Where trained networks are stored
* `notebooks/`: Collection of Jupyter Notebooks for experimenting with approximating droplet
  characteristics (data generation, training, evaluation, etc)
* `python/`: Python packages used by Jupyter Notebooks and command line scripts

# Python Setup

The `droplet_approximation` Python package can be installed into any existing
environment via `pip` though the instructions below demonstrate how to install
a Miniconda Python, create an evironment, and then install the package.  This
approach provides a consistent working environment across collaborators who may
not be regular Python developers.

* Download and Install Miniconda
* Create an NTLP environment
* Install `droplet_approximation`

## Download and Install Miniconda
Download and install Miniconda into your home directory.

***NOTE***: Miniconda is used to minimize the initial download and installation
footprint.  Anaconda may be used instead if a larger local repository of Python
packages is desired.  Replace all instances of "miniconda" with "anaconda"
below, while respecting capitalization, to switch.

[Download the appropriate client](https://www.anaconda.com/download/success).
This downloads the latest version of Python available for the x86-architecture
(since `cswarmfe` is x86-64).

```shell
$ curl -L -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Perform an unattended install into `${HOME}/miniconda3`.  Can take some time
depending on the file system's performance.

```shell
$ sh ./Miniconda3-latest-Linux-x86_64.sh -b
```

Update the shell environment to load Miniconda by default.  This updates
`${HOME}/.bashrc`.

***NOTE***: Ignore the instructions to close and re-open your current shell.
We'll do that later after additional configuration.

```shell
$ ${HOME}/miniconda3/bin/conda init
```

Configure Miniconda to use both the default and the Conda Forge channels to
locate and install packages.

```shell
$ cat > ~/.condarc <<EOF
channels:
  - conda-forge
  - defaults

channel_priority: true
EOF
```

## Create an NTLP Environment
In general, do not install new packages into the default (`base`) environment,
but rather clone it and add packages there.  We create a `ntlp` environment
and use that until a new environment is needed.

```shell
# NOTE: This can take a while!
$ conda create -n ntlp --clone base
```

Update the shell configuration so the `ntlp` environment is used.

```shell
$ echo "conda activate ntlp" >> ${HOME}/.bashrc
```

Close the current shell and re-open it (i.e. log out and log back in).  You
should now see `(ntlp)` at the beginning of your shell prompt like so:

``` shell
(ntlp) [user@system ~]$
```

## Install `droplet_approximation`'s Dependencies
The `droplet_approximation` package's dependencies are installed in two parts:

- PyTorch dependencies
- Everything else

First we install everything except PyTorch (this can take 10-15 minutes):

```shell
$ conda install --file ${HOME}/code/NTLP/droplet_approximation/python/conda_packages.txt
```

While we use Miniconda to manage our Python environments, dealing with PyTorch's
dependencies is a headache as it no longer supports installation via Miniconda.
Additionally, care should be taken to correctly choose whether one needs PyTorch
to be GPU accelerated or if it should only utilize the system's CPU(s).  While
the GPU version of PyTorch can fallback to CPU execution this is both
inefficient at run-time (i.e. not utilizing available acceleration libraries) as
well as during install-time (i.e. requiring download and installation of very
large GPU-only packages).

If CPU-only execution is desired, install PyTorch like so:

```shell
$ pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

Otherwise, install a GPU-enabled version of PyTorch:

```shell
$ pip3 install torch torchvision torchaudio
```

## Install the `droplet_approximation` Package
While we install the `droplet_approximation` package using `pip` we install its
dependencies using Miniconda to reduce headaches induced by mixing the two.

While the `droplet_approximation`'s dependencies are mostly managed by
Miniconda, the package itself is installed via `pip`.  The following
installs a version that can be modified after it has been installed so that
changes made in the NLTP working copy can be executed without reinstalling
the package.

```shell
$ pip install -e ${HOME}/code/NTLP/droplet_approximation/python
```

***NOTE***: One must reload the `droplet_approximation` module via the
`importlib` package when developing in a Jupyter Notebook as Python packages are
cached once they are initially imported.  The following code can be executed
to force a reload:

```python
import importlib
importlib.reload( droplet_approximation )
```
