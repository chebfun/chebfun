function out = feval(f, x, varargin)
%FEVAL   Evaluate the specified function.
%   Y = FEVAL(F, X) evaluates the UNBNDFUN F at the points X. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Map the input:
z = f.mapping.Inv(x);

% Make sure -Inf and Inf are mapped to -1 and 1 respectively:
mask = isinf(x); 
z(mask) = sign(x(mask));

% Evaluate the ONEFUN:
out = feval(f.onefun, z);
    
end
