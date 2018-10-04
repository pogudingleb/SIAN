# Identifiability

Maple code for assessing identifiability/observability (local and global) for models defined by systems of ODEs presented. Mostly based on the paper [Global Identifiability of Differential Models](https://cs.nyu.edu/~pogudin/global.pdf).
Tested with Maple 2018, Maple 2017, and Maple 2016.

## Usage
The main function is **IdentifiabilityODE(system, parameters, p, infolevel, method, num_nodes)** with arguments
 * **system** - a system of ODEs in the state-space form. It should include equations of two types
   * *ODEs* with rational right-hand side defining the evolution of the state variables
   * equations of the form *output_variable = rational_function(state_variables, parameters, inputs)* defining the output variables
 * **parameters** - a list of parameters and initial values whose identifiability is to be assessed. You can use **GetParameters(system)** if you want to check identifibaility of all the parameters and initial values
 * **p** - the probability of correctness, the default value is 0.99
 * **method** - the method of checking the consistency in Step 4 of Algorithm 1 from the paper. Possible options are
   * **1** - using saturation and Groebner bases, see item (1) in Remark 7 from the paper, this is usually faster
   * **2** - without saturation, with checking memebership using Groebner bases, see item (2) in Remark 7 from the paper
 * **infolevel** - the variable that regulates the amount of information printed. Options are (the default value is **1**)
   * **0** - nothing is printed
   * **1** - information about the original system and the summary of the results are printed
   * **2** - debugging mode, a lot of information is printed

Examples of usage can be found in the **examples/** folder

## Files

* **GlobalIdentifiabiliyODE.mpl**   contains the algorithm
* **/examples**   folder contains examples
  
The software is partially supported by the National Science Foundation.
