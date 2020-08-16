# finite-element-based-structural-non-linear-analysis

## Introduction ##
This is an attempt to implement a CLI tool for the following paper.

*M.C.M. Rajapakse, K.K. Wijesundara, R. Nascimbene, C.S. Bandara, R. Dissanayake, **Accounting axial-moment-shear interaction for force-based fiber modeling of RC frames, Engineering Structures**, Volume 184, 2019, Pages 15-36,ISSN 0141-0296,
https://doi.org/10.1016/j.engstruct.2019.01.075.*


**People:** [Pubudu Premathilaka], [Suneth Samarasinghe]

**Advised by:** [Dr. Kushan Wijesundara](http://eng.pdn.ac.lk/civil/people/drKKWijesundara.php) and [Sameera Hippola)


## Project summary ##

* Civil engineering structure modelling
  * Elementwise local stiffness matrix calculation generation.
  * Structure's global stiffnexx matrix calculation.
  * Load matrix generation.
* Lieanr system oprations.
  * Linear system representation : dense and sparse
  * Linear system solving: Gauss ellimination on dense and sparse matrices, iterative numerical techniques.
* Backtracking to translate the linear system solution to civil engineering structure.

The main objective of this project was to find the **time and memroy efficient** way of computing the results.


## Please note ##
This work was done as a partial requirement for [CO328 Software Engineering](http://www.ce.pdn.ac.lk/undergraduate-courses/)) course. **This implementation consists of the linear and non linear region analysis **

