function C = plus(A, B)
%+   CHEBOP plus.
%   C = A + B, where A and B are CHEBOP objects, returns a CHEBOP C
%   representing the addition of the operators of A and B. Boundary
%   conditions on A or B are maintained, but only one of A and B should
%   have boundary conditions. The .op fields of A and B must have the same
%   input variables.
% 
% Example:
%
%   L = chebop(@(u) -diff(u,2), [0,pi], 'dirichlet');
%   I = chebop(@(u) u, [0,pi]);
%   eigs(L+I)
%
% See also CHEBOP/MINUS, CHEBOP/MTIMES, CHEBOP/MLDIVIDE, CHEBOP/FEVAL

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(A,'chebop') && isa(B,'chebop') )       
    
    if ( hasbc(A) && hasbc(B) )
        error('CHEBFUN:CHEBOP:plus:bcs', ...
        'Only one input may have boundary conditions in CHEBOP + CHEBOP.');
    end
    
    % Initialise output (make sure it contains the BCs)
    if ( hasbc(A) )
        C = A;
        A = B;
    else
        C = B;
    end

    % Create new anon. func
    funArgs = getFunArgs(C);
    if ( ~strcmp(funArgs, getFunArgs(A)) )
        error('CHEBFUN:CHEBOP:plus:inputArgs', ...
        'Incompatable input arguments in A.op and B.op.');
    end
    C.op = eval(['@(',funArgs,') C.op(',funArgs,') + A.op(',funArgs,')']);         
    
else 
    error('CHEBFUN:CHEBOP:plus:notSupported', ...
        '%s + %s addition is not supported.', class(A), class(B));
    
end
    
end
