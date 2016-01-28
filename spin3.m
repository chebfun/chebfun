function [uout, tout] = spin3(varargin)
%SPIN3  Solve a time-dependent PDE in 3D with periodicity in space, using a 
%Fourier spectral method and an exponential integrator time-stepping scheme.
%   SPIN3

% See also SPINOP3, SPINPREF3, SPINSCHEME, SPIN, SPIN2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN3 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end