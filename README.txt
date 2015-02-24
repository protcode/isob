Thank you for downloading the isobarQuant package. By using this package you agree to abide by the terms of the License accompanying this release. 

In order for it to run you will need to install the Python programming language (for windows) and some extra Python packages as described below.

Please note that this package is designed to run on the Windows operating system (it is possible to run the postMascot workflow on Linux but it has not been tested on other operating systems). The preMascot section of the workflow also requires the Xcalibur Foundation software to be installed and has been tested with Xcalibur versions 2.2 and 3.0. 

The package has been developed using Python2.7.1 but will run with any Python2.7 version. It is not yet compatible with Python3x. Note that Xcalibur is limited to a 32bit Windows environment: it was only possible to compile the Xcalibur/Python API to run with a Python 32bit version. 

Once you have completed the steps below you should be ready to run the package. If the machine on which the installation is taking place has an http(s) connection to the outside world follow steps 1-3, if you are behind a firewall or have no access to download packages, steps 1-3 should be replaced by i to iv, where the packages are downloaded on another computer and then copied to a local directory (eg c:\temp) from which Python can install them. If using an earlier version of Python than 2.7.9 it will be necessary to install the Python 'pip' package as well. Ways to do this can be found here: https://pip.pypa.io/en/latest/installing.html 

The example below is for 2.7.9 and assumes that c:\temp is the location of all files downloaded. 

1) download the following Microsoft installer (msi) file:  https://www.python.org/ftp/python/2.7.1/python-2.7.9.msi 
Execute the .msi file and follow instructions to install Python to a location of your choice. Ensure that Python is added to the system path (this is the bottom option on the first page of the installer after the point where the user is requested to set the location for python to be installed)

2) From http://www.lfd.uci.edu/~gohlke/pythonlibs/ download the following 3 binary (.whl) files, which correspond to Python version 2.7 on a 32-bit Windows machine to any (temporary) folder eg c:\temp\ 

numpy-1.8.2+mkl-cp27-none-win32.whl 
numexpr-2.4-cp27-none-win32.whl 
tables-3.1.1-cp27-none-win32.whl

and install them via the command prompt in the order given below:

 C:\>pip install c:\temp\numpy-1.8.2+mkl-cp27-none-win32.whl 
 C:\>pip install c:\temp\numexpr-2.4-cp27-none-win32.whl 
 C:\>pip install c:\temp\tables-3.1.1-cp27-none-win32.whl

3) From a command prompt run: C:\pip install pathlib

Installation finished - the downloaded files in C:\temp may now be deleted

#--------------------------------------------#

Alternatively, if the computer onto which you are installing the packages is behind a firewall / unable to connect to the internet the following procedure should be followed:

i) On a computer with internet connection, download the following Microsoft installer (msi) file:  https://www.python.org/ftp/python/2.7.1/python-2.7.9.msi and copy it to a directory of your choosing on the computer where you are installing isobarQuant. eg c:\temp

ii) On a computer with internet connection, download the following 3 binary (.whl) files, from http://www.lfd.uci.edu/~gohlke/pythonlibs/:
numpy-1.8.2+mkl-cp27-none-win32.whl 
numexpr-2.4-cp27-none-win32.whl 
tables-3.1.1-cp27-none-win32.whl
and also copy them to directory of your choosing on the computer where you are performing the installation eg c:\temp

ii) On a computer with internet connection, download pathlib from:
https://pypi.python.org/packages/source/p/pathlib/pathlib-1.0.1.tar.gz
and copy to a directory of your choosing on the computer where the installation is being done eg c:\temp 

iii)Execute the .msi file and follow instructions to install Python to a location of your choice. Ensure that Python is added to the system path (this is the bottom option on the first page of the installer after the point where the user is requested to set the location for python to be installed)

iv) On the computer where the installation is being carried out, open a command prompt and run 
  iv.1) C:\>pip install c:\temp\numpy-1.8.2+mkl-cp27-none-win32.whl
  iv.2) C:\>pip install c:\temp\numexpr-2.4-cp27-none-win32.whl
  iv.3) C:\>pip install c:\temp\tables-3.1.1-cp27-none-win32.whl
  iv.4) C:\>pip install c:\temp\pathlib-1.0.1.tar.gz



Installation finished - the downloaded files in C:\temp may now be deleted

