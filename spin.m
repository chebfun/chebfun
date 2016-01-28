function [uout, tout] = spin(varargin)
%SPIN  Solve a time-dependent PDE in 1D with periodicity in space, using a 
%Fourier spectral method and an exponential integrator time-stepping scheme.
%   SPIN

% See also SPINPREF, SPINSCHEME, SPIN2, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end