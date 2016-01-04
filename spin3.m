function [uout, tout] = spin3(varargin)
%SPIN3  Solve a (t,x,y,z)-PDE with periodicity in space, using a Fourier 
%spectral method for x, y and z, and an exponential time-differencing scheme for 
%t.
%   SPIN3

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[uout, tout] = spinoperator.solvepde(varargin{:});

end