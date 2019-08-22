# SIAN v0.5

Maple implementation of the algorithm for checking global identifiability for systems of ODEs presented in the paper [Global Identifiability of Differential Models](https://cs.nyu.edu/~pogudin/global.pdf).
Tested with Maple 2017 and Maple 2016.

## Usage
The main function, which corresponds to Algorithm 1 from the paper, is **GlobalIdentifiability(sigma, theta_l, p, method)** with arguments
 * **sigma** - a Maple table that describes a system of ODEs with the following entries
   * **x_vars** - the list of state variables (**x** in the paper)
   * **x_eqs** - the list of state equations (**x' = f(x, mu, u)** in the paper) in the *same* order as the variables in x_vars
   * **y_vars** - the list of output variables (**y** in the paper)
   * **y_eqs** - the list of output equations (**y = g(x, mu, u)** in the paper) in the *same* order as the variables in y_vars
   * **u_vars** - the list of input variables (u in the paper)
   * **mu** - the list of parameters that are *not* the initial conditions (**mu** in the paper)
 * **theta_l** - a subset of locally identifiable parameters
 * **p** - the probability of correctness, the default value is 0.99
 * **method** - the method of checking the consistency in Step 4 of Algorithm 1 from the paper. Possible options are
   * **1** - using saturation and Groebner bases, see item (1) in Remark 7 from the paper, this is usually faster
   * **2** - without saturation, with checking memebership using Groebner bases, see item (2) in Remark 7 from the paper

Examples of usage can be found in the **examples/** folder

## Files

* **GlobalIdentifiabiliy.mpl**   contains the algorithm
* **/examples**   folder contains examples from the paper
  * **example_paper.mpl** contains the code for Example 5
  * **ChemicalReactionNetwork.mpl** contains the code for Example 6
  * **Cholera.mpl** contains the code for Example 7
