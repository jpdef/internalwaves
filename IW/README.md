# Internal Waves
Library of routines to calculate internal wave fields, and invert internal wave field data make to model paramters

## Running Simulations
The simulations take a configuration file which specifies the paramters of the simulation. They are written in JSON format. Examples are found in ./script/config. 
To execute a simulation run the following command. 
```
cd scripts
./simulate.py <config file>
```
This wil output either data files ( feather files) or a movie to the specified path under the data folder.

##Data Output
The data files are binary files that represent a time frame. Each frame has a set of rows which contain the unique position , time , and displacement value. The file size will scale by the amount of spacial positions the user choices to sample (num x) * (num y ) * (num z). There is also a meta data file stored in the data directory that contains various parameters of the simulation. Each feather file can be read in as a dataframe using each python (panda dataframes) or R native dataframes.

```
import pandas 
df = pandas.read_feather('data/myexample.fthr')
```


