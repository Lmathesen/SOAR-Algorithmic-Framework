# SOAR-Algorithmic-Framework
Stochastic Optimization with Adaptive Restart (SOAR) a Global Optimization Framework

This project contains a selection of implementations of the SOAR framework. If using this code, please refer and cite the paper at the close of this ReadMe for algorithmic theoretical presentations.

To execute the code in this repository, all files and folders must necissarily be downloaded. The executable main method file "SOAR_algorithm" contains locations for algorithmic user input parameterization, as well as several sections (to be detailed below) which can be altered to deliver alternative SOAR implementations. 
The current code in the reposity is configured to match the implementation used to gather the results presented in the paper  

## Authors: 
 * **Logan Mathesen** lmathese@asu.edu - *Arizona State Univeristy*
 * **Giulia Pedrielli** giulia.pedrieli@asu.edu - *Arizona State Univeristy* 
     
## Algorithm and Code Interface     
As outlined in the paper, this algorithm contains three core components.

### Global Surrogate Modeling
	* We make use of spatial Gaussian Process Models, specifically the Ordinary Kriging form with diagonal covariance matrix perturbations for numerical parameter estimation.
	* All files needed for global surrogate estimation are contained in the "Ordinary Kriging Implementation" folder.
	* All coded algorithms in this repository execute a constant mean form with Gaussian Correlation function, please refer to the files in this folder for alternative Kriging model implementations for different modeling assumptions.

### Local Search Component
* We make use of a Trust Region Subproblem based approach for local optimization. In the main method 'SOAR_algorithm' file, this entire procedure (with adaptive restart identification integrated) composes lines 138 through 355.  
* There are two alternatives for the trust region subproblem solved, the quadratic and linear taylor expansion local model forms about the current centroid.
** Both the linear and quadratic model form function files reside in the "Trust Region Implementation" folder, along with the gradient and hessian estimations function files.  
** Currently the quadratic model TR subproblem is set as the default alternative, this can be switched to the linear subproblem by altering lines 276 and 280. 
*** Lines 276 and 280, in three locations total, currently read ...quadratic_model(.,.,G1',H)...
*** All three of these locations should be altered to read ...linear_model(.,.,G1')... to enable the linear model, which no longer needs hessian estimates
note if this subproblem is utlized, lines 261 and 354 should be commented. Additionally the for loop condition in lines 220, 250, 318, and 343 need to be altered to "for i = 1:(problem.dim)"  
	
### Adaptive Restart Process
	* Sampling Criteria - The "Improvement Functions" folder contains two basic alternative sampling criteria function files. The first delivers expected improvement (EI) function value predictions, the second delivers the pre-integrated form of this function- which we term probabilistic improvement (PI) fucntion values.  
	* In general there are three operating modes for the code with respect to the implemented sampling criteria, follow the instructions below to alter the main method file "SOAR_algorithm" to deliver your desired adaptive restart sampling criteria implementation.
		** Crowded-Expected Improvement: The code in the repository is set to this default configuration. The function in line 116 must read EI_calc_kd_pred(...). Lines 120-131 and 135-148 must be commented out. Section in lines 153-182 select the next restart centroid as that which maximizes crowding distance subject to belong to the top alpha-level set (this section to remain commented as is). Lines 185, 186 must be commented.
		** Expected Improvement: The function in line 116 must read EI_calc_kd_pred(...).  Lines 120-131 and 135-148 must be commented out. Entire section, lines 153 - 182 must be commented out. Lines 185 and 186 must be uncommented.
		** Probabilistic Improvement: The function in line 116 must read [PI_values, Varextrinsic, y_fit] = PI_calc_kd_pred(...).  Lines 120-131 and 135-148 must be uncommented. Entire section, lines 153 - 182 must be commented out. Lines 185 and 186 must be commented.
	* Dynamic Restart criteria is coded as presented in the paper, and can be found in lines 288-292. There are currently no alternatives coded for this criteria. 
		
	
    
## User defined inputs for overall SOAR algorithm execution (lines 14-47).
	* Altering test functions environment. 
		** Uncomment the line of the test problem you wish to solve, comment all the others. All preloaded test functions are stored in the "Test Cases" folder. 
		** Follow/Drill into the problem class for your desired problem to alter problem properties, properties include: problem dimension "object.dim," and upper/lower bounds "obj.ub" and "obj.lb." Note only the kd (k-dimensional) problems allow specificiation of any problem dimension.  
	* Altering Discritized Design Grid.
		** Note the discritized design grid must match dimension with the specified problem to be solved. 
		** 3 designs grids (4, 8, and 17 dims) are preconstructed, to reproduce the results derived in the paper. These can be enabled by commenting lines 45,46 and uncommenting line 44 substiting the appropriate file to be loaded 'design_grid(100000x17)' for 17 dims, 'design_grid(75000x8)' for 8 dims, and 'design_grid(40000x4)' for 4 dims. These .mat data files can be found in the "Cross Validation and Initialization" folder in the repository.
		** If choosing your own dimension, set the size of the design grid appropriately (line 45). In line 46, the 'Cdf' function (in folder "Cross Validation and Initialization") will construct the corresponding grid for your problem. 
	* Crowded EI level set threshold can be adjusted, line 32, and Trust Region control and testing parameters can be adjusted in lines 35-40.

## Execution Result Storage
Note that at the close of the main method the results of the algorithmic experimentation can be exported via the 'x_storage' variable (where the minimum observed function value location is stored immediately following each function value observation). This variable stores algorithm replications in rows, each simulation record in columns, and each dimension of that given record in pages. 
	If you would like to export your resutls for further processing and analysis, please enter your file export location in the xlswrite3 function.


## Please cite the following papers if you use this code: 
```
Mathesen, L., Pedrielli, G. (2018). Stochastic Optimization With Adaptive Restart: A Framework for Integrated Local and Global Learning. Manuscript submitted for publication.
```
```
Mathesen, L., Pedrielli, G., & Ng, S. H. (2017, December). Trust region based stochastic optimization with adaptive restart: A family of global optimization algorithms. In Simulation Conference (WSC), 2017 Winter (pp. 2104-2115). IEEE.
```
