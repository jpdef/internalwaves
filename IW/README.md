# internalwaves
Library of routines to calculate internal wave fields, and invert internal wave field data make to model paramters

##Running Simulations
The simulations take a configuration file which specifies the paramters of the simulation. They are written in JSON format. Examples are found in ./script/config. 
To execute a simulation run the following command. 
```
./run_sim.py <config file>
```
This wil output either data files ( feather files) or a movie to the specified data path in the configuration file.


