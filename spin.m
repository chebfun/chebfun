function [uout, tout] = spin(varargin)
%SPIN  Solve a time-dependent PDE in 1D with periodicity in space, using a 
%Fourier spectral method in space and an exponential integrator in time.
%   SPIN

% See also SPINPREF, SPIN2, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end