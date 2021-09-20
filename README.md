# boolean_circuit_model

C++ source code of the Boolean circuit model used in:

Wechsler, D., and J. Bascompte. 2021. Cheating in mutualisms promotes diversity
and complexity. The American Naturalist, in press.

Some of the parameters in the program code are not documented in the manuscript. For the simulation experiments done in the manuscript the default values of these parameters were used. 

Further, there are differences with regard to the names of certain model parameters and variables between the model description in the manuscript and the code. Specifically, these differences are:

* Within the program animals are called active agents and plants passive 
  agents.

* The variable `n` (number of traits/inputs) in the code is incremented by 
  2 compared to the manuscript (because in the code the 
  two constant inputs also count as inputs). Hence, to simulate a 
  scenario with 3 traits use `n=5` in the code.

* The number of rows of the circuit template is called `r` in the code 
  (instead of `k`)

* The number of columns of the circuit template is called `c` in the code
  (instead of `m`)

* The epistatic function is called `f_passive` in the code (instead of `f`).


## REQUIREMENTS

In order to build `cheater_recognition`, you must have [CMake](http://www.cmake.org) installed.

Further, the program uses the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/).

## INSTALLATION AND EXECUTION

To compile and install the program, run the following commands 
in a terminal (in the root directory of the progject):

`$ cmake .`

`$ make`

`$ make install`

A new directory called 'bin' is created. It contains the executable 
`cheater_recognition` which can be started using the following command
(from within the `bin` directory):

`$ ./cheater_recognition`

This will start the program using default parameters. To find out about 
the possible parameters run the command: 

`$ ./cheater_recognition --help`

(See also Table S1 in the supplementary material):


## SIMULATION OUTPUT

Simulation results are written to the standard output in CSV format. It contains 
for each `SAMPLE_INTERVAL`-(th) generation between `SAMPLE_START` and 
`SAMPLE_END` a row with statistics about the animal and plant population at 
that generation.

In particular, the CSV contains the following columns:

Column        | Content           
------------- |------------- 
`T` | The time step (generation).
`N_passive_X` | The number of plant individuals of type `X`. Where `X` is an integer in `{0, (2^n)-1}` that identifies the plant type (see table below).
`payoff_passive_X` | The payoff of plants of type `X`.
`N_active_Y` | The number of animals with phenotype `Y`. Where `Y` is an integer in `{0, (2^(2^n))-1}` that identifies the phenotype (see table below). There exists only a column for a phenotype if it appeared at least once during the simulation.		
`avg_payoff_active_Y` | The average payoff of animals with phenotype `Y`. The payoff of animals with the same phenotype can be different if the parameter `d` is `> 0`. The parameter `d` specifies the cost for a logic gate. For the submitted manuscript we considered only `d=0`.
`avg_num_gates_Y` | The average number of logic gates of the animals with phenotype `Y`.

The following table illustrates, by means of an example with 3 traits, 
how the integer numbers used to identify plant types (`X`) and animal 
phenotypes (`Y`) in the output CSV are encoded.

The columns `t1`, `t2` and `t3` list all possible plant types for 3 traits.
The traits are interpreted as a binary number with `t3` the least and `t1` 
the most significant bit. The integer representation of the plant type in 
the left most column is the decimal conversion of this binary number.

The right column lists some of the possible phenotypes (0=avoid, 1=interact).
The phenotype is also interpreted as a binary number with the first row 
holding the least and the last row the most significant bits. The integer 
representation of the phenotype at the bottom is the decimal conversion of 
this binary number.

	+-----------+-----------------+-------------------------+
	|  plant    |   plant traits  |     animal phenotypes   |
	|  decimal  |    t1  t2  t3   |                         |
	+-----------+-----------------+-------------------------+
	|    0	    |    0   0   0    |   0  1  0  ..  1 .. 1   |
	|    1	    |    0   0   1    |   0  0  1  ..  1 .. 1   |
	|    2	    |    0   1   1    |   0  0  0  ..  1 .. 1   |
	|    3	    |    0   1   1    |   0  0  0  ..  1 .. 1   |
	|    4	    |    1   0   0    |   0  0  0  ..  0 .. 1   |
	|    5	    |    1   0   1    |   0  0  0  ..  0 .. 1   |
	|    6	    |    1   1   0    |   0  0  0  ..  0 .. 1   |
	|    7	    |    1   1   1    |   0  0  0  ..  0 .. 1   |
	+-----------+-----------------+-------------------------+
	|  	       animal decimal ->  0  1  2      15   255 | 
	+-------------------------------------------------------+









