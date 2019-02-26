## !!This repository and the documentation are under development)!!

# Test-particle  
This code trace the trayectory of energetic particles in the Radiation belts. The method used is described in 
##### Fok, M.-C., and T. E. Moore. 1997. "Ring current modeling in a realistic magnetic field configuration." Geophys Res Lett, 24: 1775-1778 [10.1029/97GL01255]

The original test-particle code was written in fortran by Fok, M.-C. and collaborators. This project convert the original code to python and improve the visualitation and analysis of the data. 

Also it is added the possibility of tracing the trayectory of particles with non 90 degrees pitch angles.
## Dependencies

AE index plot (https://github.com/scivision/AEindex#ae-index-plot)

spacepy https://pythonhosted.org/SpacePy/install.html

## Install

1- Clone and intall

git clone https://github.com/germanfarinas/test-particle.git
cd test-particle
python setup.py install


2- Dowload the the Qin-Denton database

The Tsyganenko routines use the W parameters. These parameters are already pre-calculated in the Qin-Denton database. To download the data execute the following command


3- 
