function L = lagpoly( n )
% LAGPOLY   Laguerre polynomial of degree n.
% L = LAGPOLY(N) returns the chebfun corresponding to the Laguerre polynomials 
% L_N(x) on [0,inf], where N may be a vector of positive integers.
%
% Note, this is currently just a toy to play with the construction of
% Hermite polynomials using a combination of Chebfun's barycentric,
% mapping, and 'blowup' technologies. See chebfun/chebtests/unbndpolys.m
% for some testing.
%
% See also chebpoly, legpoly, jacpoly, and hermpoly.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

xinf = chebfun(@(x) x,[0 inf]);
L = chebfun(@(x) 1+0*x,[0 inf]);
L = [L chebfun(@(x) 1-x, [0 inf] )];
for k = 2:max(n) % Recurrence relation
    Lk = ((2*k-1-xinf).*L(:,k) - (k-1)*L(:,k-1))/k;
    L = [L Lk]; 
end

% Take only the ones we want
L = L(:,n+1);
