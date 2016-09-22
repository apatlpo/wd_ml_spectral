
This library solves the response a mixe layer to a time periodic wind on
a spatially uniform horizontal grid
It relies on petsc4py for parallelized inversion 

Download and install with:
```
git clone https://apatlpo@bitbucket.org/apatlpo/wd_solver.git
cd wd_solver 
python setup.py build_ext --inplace 
```


We use conda for the python install:
```
bash
source activate natl60
```

Proper conda install on Linux:
```
bash /home/mulroy/slgentil/tarfiles/Miniconda2-latest-Linux-x86_64.sh
bash
conda update conda
conda create --name natl60 python
source activate natl60
conda install dask
conda install xarray
conda install -c juanlu001 petsc4py=3.6.0
conda install libgfortran=1.0
conda install -c scitools cartopy=0.13.1
conda install basemap
conda install -c asmeurer pango
```

Proper conda install on Caparmor:
```
conda create --name petsc python
source activate petsc
conda install -c sed-pro-inria petsc4py=3.4
conda install -y netcdf4
```
