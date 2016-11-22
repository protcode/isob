# isobarquant

Thank you for downloading the isobarQuant package. By using this package you agree to abide by the terms of the License accompanying this release.

## Installation

In order for it to run you will need to install the Python programming language (for windows) and some extra Python packages as described below.

Please note that this package is designed to run on the Windows operating system (it is possible to run the postMascot workflow on Linux but it has not been tested on other operating systems). The preMascot section of the workflow also requires the Xcalibur Foundation software to be installed and has been tested with Xcalibur versions 2.2 and 3.0.

The package has been developed using Python2.7.1 but will run with any Python2.7 version. It is not yet compatible with Python3x. Note that Xcalibur is limited to a 32bit Windows environment: it was only possible to compile the Xcalibur/Python API to run with a Python 32bit version.

(The example below is for 2.7.9 - if using an earlier version of Python than 2.7.9 it will be necessary to install the Python 'pip' package as well. Ways to do this can be found here: https://pip.pypa.io/en/latest/installing.html)

We recommend removing any previous installations of Python before installing the new version detailed below.

1) Download the following Microsoft installer (msi) file:  https://www.python.org/ftp/python/2.7.9/python-2.7.9.msi
Execute the .msi file and follow instructions to install Python to a location of your choice. Ensure that the location of the Python interpreter (python.exe, probably located in C:\Python27\) is added to the system path (this is performed during installation by selecting the bottom option on the first page of the installer after the point where the user is requested to set the location for python to be installed, by default it is set as a red 'X')

2) Within this download there are 3 binary (.whl) files taken from the external repository maintained by Christoph Gohlke (http://www.lfd.uci.edu/~gohlke/pythonlibs/), which correspond to Python version 2.7 on a 32-bit Windows machine:

- numpy-1.8.2+mkl-cp27-none-win32.whl
- numexpr-2.4-cp27-none-win32.whl
- tables-3.1.1-cp27-none-win32.whl

3) Use `pip` to install these packages/files like so: 

    pip install numpy-1.8.2+mkl-cp27-none-win32.whl
    pip install numexpr-2.4-cp27-none-win32.whl
    pip install tables-3.1.1-cp27-none-win32.whl

4) Finally install the pathlib library either from the web:

    pip install pathlib==1.0.1

  or by using the provided archive:

    pip install pathlib-1.0.1.tar.gz
