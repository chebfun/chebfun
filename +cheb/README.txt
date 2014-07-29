The purpose of this package is to provide extra capabilities not found in the
core Chebfun classes and functions. To be added to this folder, functions 
must have the following properties:

1. No other functions can depend on it.

2. It depends only on built-in MATLAB and core Chebfun. 

3. It is well documented, following MATLAB conventions for the help 
documentation block. In particular, it should provide at least one example 
of usage that will be printed out by “help”.

4. It has a comment with the author’s name and email, and the month and year 
of original authorship. 

5. It has a notice assigning copyright to The University of Oxford and The 
Chebfun Developers. However, you should include your own 

Submit a pull request for your function, adding a justification for the 
existence of the function. The code will be subject to review by at least 
one developer. 