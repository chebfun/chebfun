function A = uminus(A)
%-   CHEBOP unary minus.
%   -A negates the operator of the CHEBOP A. Boundary conditions are not
%   affected.
%
%   UMINUS(A) is called for the syntax '-A'.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(A.op ) )
    return
else
    % Multiplication of anonymous functions is not supported in Matlab. 
    % Need to work around that.
    funArgs = getFunArgs(A);
    A.op = eval(['@(', funArgs, ') -A.op(', funArgs, ')']);
end
    
end
