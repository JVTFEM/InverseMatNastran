# About
This repository contains all the source code and an example problem to perform material characterization from digitial image correlation (DIC) data.  The process makes use of Nastran to perform the FE analyses and DOT (Design Optimization Tools from [VR&D](http://www.vrand.com)) as optimizer.  It is a gradient-based optimizer.

# Information
This project was setup and tested on a Linux environment, but should also run an Windows with minor changes to the code.

## Requirements
- The DOT optimizer (commericial code form [VR&D](http://www.vrand.com)).  A Python wrapper for DOT is provided that will allow the user to call the DOT shared library from Python.
- [pyNastran](https://pynastran-git.readthedocs.io/en/latest/) for working with Nastran input/output files.
- [Pandas](http://pandas.pydata.org/) Python library
- [NumPy](http://www.numpy.org/) Python library

# How it works
The process required to evaluate the objective funcstion can be summarrize das:

# File description


# How to run
Before running the code, please make sure that you haveDOT instlaled with a license file.  Alsoe make sure that you have Nastran installed and know twhwere the path is is.

1. Clone this library
'git clone https://github.com/cjekel/shape_optimization_rubber_gasket_Abaqus.git'

2. Move into the shape_opt folder
'cd shape_opt'

3. Chagne the path for NASTRAN and gzip

3. Run the optimization
'python optimizatin.py'
