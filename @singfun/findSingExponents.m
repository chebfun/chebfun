function exponents = findSingExponents(f, singType)
%FINDSINGEXPONENTS   Endpoint singularity detection by sampling values.
%  Private method of SINGFUN.

% TODO Documentation

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

tol = singfun.pref.singfun.eps;

x1 = eps; x2 = 2*x1;
exponents(1) = log(f(-1+x2)/f(-1+x1))/log(2);
exponents(2) = log(f(1-x2)/f(1-x1))/log(2);