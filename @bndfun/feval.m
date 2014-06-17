function out = feval(f, x, varargin)
%FEVAL   Evaluate the specified function
%   Y = FEVAL(G, X) Evaluation of a BNDFUN G at points X. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Map the input:
z = f.mapping.Inv(x);

% Evaluate the onefun:
out = feval(f.onefun, z);
    
end
