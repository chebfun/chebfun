function C = plus(A, B, p)
%+   CHEBOP plus.
%   C = A + B, where A and B are CHEBOP objects, returns a CHEBOP C
%   representing the addition of the operators of A and B. Boundary
%   conditions on A or B are maintained, but only one of A and B should
%   have boundary conditions. The .op fields of A and B should have the
%   same input variables, but this is not validated with this PLUS method.s
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
        flip = true;
    else
        C = B;
        flip = false;
    end

    funArgs1 = getFunArgs(C);
    funArgs2 = getFunArgs(A);
    argStr1 = ['@(',funArgs1,')'];
    argStr2 = ['@(',funArgs2,')'];
    
    % Pretty print:
    if ( ~isempty(C.opShow) )
        op = C.opShow;
    else
        op = func2str(C.op);
    end
    op = strrep(op, argStr1, '');
    if ( ~isempty(A.opShow) )
        op2 = A.opShow;
    else
        op2 = func2str(A.op);
    end
    op2 = strrep(op2, argStr2, '');
    if ( nargin == 2 )
        p = '+';
    end
    if ( flip )
        C.opShow = ['@(',funArgs1,')' op p op2];
    else
        C.opShow = ['@(',funArgs1,')' op2 p op];
    end
    
    % Create new anon. func:
    C.op = eval(['@(',funArgs1,') C.op(',funArgs1,') + A.op(',funArgs1,')']);         
    
else 
    error('CHEBFUN:CHEBOP:plus:notSupported', ...
        '%s + %s addition is not supported.', class(A), class(B));
    
end
    
end
