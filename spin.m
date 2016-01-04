function [uout, tout] = spin(varargin)
%SPIN  Solve a (t,x)-PDE with periodicity in space, using a Fourier spectral 
%method for x and an exponential time-differencing scheme for t.
%   SPIN

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[uout, tout] = spinoperator.solvepde(varargin{:});

end