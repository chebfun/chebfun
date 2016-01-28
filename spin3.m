function [uout, tout] = spin3(varargin)
%SPIN3  Solve a time-dependent PDE in 3D with periodicity in space, using a 
%Fourier spectral method in space and an exponential integrator in time.
%   SPIN3

% See also SPINPREF3, SPIN, SPIN2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN3 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end