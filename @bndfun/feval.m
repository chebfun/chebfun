function out = feval(f, x)
%FEVAL   Evaluate the specified function
%   Y = FEVAL(G, X) Evaluation of a BNDFUN G at points X. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Map the input:
z = f.mapping.inv(x);

% Evaluate the onefun:
out = feval(f.onefun, z);
    
end
