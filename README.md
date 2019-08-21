# SIAN (Structural Identifiability ANalyser)

Maple code for assessing identifiability/observability (local and global) for models defined by systems of ODEs presented. Mostly based on the paper [Global Identifiability of Differential Models](https://cs.nyu.edu/~pogudin/global.pdf). Supplementary Maple code for the paper is available at https://github.com/pogudingleb/Global_Identifiability
Tested with Maple 2018, Maple 2017, and Maple 2016.

## How to download
* Usign *git checkout*
* On https://github.com/pogudingleb/SIAN, press green button "Clone of download", and then "Dowload ZIP"

## Usage
The main function is **IdentifiabilityODE(system, parameters, p, infolevel, method, num_nodes)** with required positional arguments **system** and **parameters** and optional keyword arguments **p**, **infolevel**, **method**, **num_nodes**.
 * **system** - a system of ODEs in the state-space form. It should include equations of two types
   * *ODEs* with rational right-hand side defining the evolution of the state variables
   * equations of the form *output_variable = rational_function(state_variables, parameters, inputs)* defining the output variables
 * **parameters** - a list of parameters and initial values whose identifiability is to be assessed. You can use **GetParameters(system)** if you want to check identifibaility of all the parameters and initial values
 * **p** (optional) - the probability of correctness, the default value is 0.99. For technical reasons
 * **infolevel** (optional) - the variable that regulates the amount of information printed. Options are (the default value is **1**)
   * **0** - nothing is printed
   * **1** - information about the original system, main steps of the algorithm, and the summary of the results are printed
   * **2** - debugging mode, a lot of information is printed
 * **method** (optional) - the method of checking the consistency in Step 4 of Algorithm 1 from the paper. Possible options are (the default value is **1**)
   * **1** - using saturation and Groebner bases, see item (1) in Remark 7 from the paper, this is usually faster
   * **2** - without saturation, with checking memebership using Groebner bases, see item (2) in Remark 7 from the paper
 * **num_nodes** (optional) - the maximal number of processes created by the algorithm, the default value is 6.

*Example* **IdentifiabilityODE(s, [a, b], infolevel = 2, num_nodes = 5)**. For more details on positional and keyword arguments in Maple, see [here](https://www.maplesoft.com/support/help/maple/view.aspx?path=parameter_classes).


Examples of usage can be found in the **examples** folder. One can run an example by
  * either opening it *as a Maple worksheet* and executing it
  * or from the *command line* by 
    * going to **examples** directory
    * calling *maple name_of_file.mpl* 

## Files

* **IdentifiabiliyODE.mpl**   contains the algorithm
* **examples**   folder contains examples
  
The software is partially supported by the National Science Foundation.
