function C = mtimes(A, B)
%*   CHEBOP composition, multiplication, or application.
%   C = A*B, where either A or B is a scalar, returns a CHEBOP C representing
%   scalar multiplication of the original operator.  Boundary conditions
%   are copied from A or B to C.
%
%   If N is a CHEBOP and U a CHEBFUN or CHEBMATRIX of dimension compatible
%   with N.op, then N*U is equivalent to FEVAL(N, U).
%
%   C = A*B, where A and B are CHEBOP objects, should return a CHEBOP C
%   representing the composition of the operators of A and B. Boundary
%   conditions on A or B are destroyed by this process. Note this is not yet
%   supported.
%
% See also CHEBOP/PLUS, CHEBOP/MINUS, CHEBOP/MLDIVIDE, CHEBOP/FEVAL

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(A, 'chebfun') )
    error('CHEBFUN:CHEBOP:mtimes:invalid', 'Operation is undefined.');
    
elseif ( isa(B, 'chebfun') || isa(B, 'chebmatrix') )    
    % Let A operate on B:
    C = feval(A, B);
    
elseif ( isnumeric(B) )
    % Switch argument to make sure A is numeric:
    C = mtimes(B, A);
    
elseif ( isnumeric(A) )
    
    if ( ~isscalar(A) )
        error('CHEBFUN:CHEBOP:mtimes:nonScalar', ...
            'CHEBOP * DOUBLE multiplication is defined only for scalars.');
    end
    
    % Initialize C as a CHEBOP
    C = B;
    
    % Make the pretty string:
    funArgs = getFunArgs(C);
    argStr = ['@(',funArgs,')'];
    if ( ~isempty(C.opShow) )
        op = C.opShow;
    else
        op = func2str(C.op);
    end
    C.opShow = strrep(op, argStr, [argStr num2str(A) '*']);
    
    % Multiplication of anonymous functions is not supported in Matlab. Need to
    % work around that.
    C.op = eval(['@(', funArgs, ') A*C.op(', funArgs, ')']);
    
elseif ( isa(A,'chebop') && isa(B,'chebop') )       
    % CHEBOP composition is not yet supported. It will probably be a mess once
    % we start getting systems involved, as the inputs need to be shuffled
    % correctly. A remedy might be through a nested function in this file.
    error('CHEBFUN:CHEBOP:mtimes:notSupported1',...
        'CHEBOP composition is not supported.');   
    
else 
    error('CHEBFUN:CHEBOP:mtimes:notSupported2', ...
        '%s * %s multiplication is not supported.', class(A), class(B));
    
end
    
end
