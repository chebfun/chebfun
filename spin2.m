function [uout, tout] = spin2(varargin)
%SPIN2  Solve a (t,x,y)-PDE with periodicity in space, using a Fourier spectral 
%method for x and y, and an exponential time-differencing scheme for t.
%   SPIN2

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[uout, tout] = spinoperator.solvepde(varargin{:});

end