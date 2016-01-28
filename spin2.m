function [uout, tout] = spin2(varargin)
%SPIN2  Solve a time-dependent PDE in 2D with periodicity in space, using a 
%Fourier spectral method in space and an exponential integrator in time.
%   SPIN2

% See also SPINPREF2, SPIN, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN2 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end