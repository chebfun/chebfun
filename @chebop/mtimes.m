function C = mtimes(A, B)
%* Chebop composition, multiplication, or application.
% If A and B are chebops, then C = A*B is a chebop where the operator of C
% is the composition of the operators of A and B. No boundary conditions
% are applied to C.
%
% If either A or B are scalar, then C = A*B is a chebop representing scalar
% multiplication of the original operator. In this case, boundary
% conditions are copied into the new operator.
%
% If N is a chebop and U a chebfun, then N*U applies N to U.
%
% See also CHEBOP/MLDIVIDE, CHEBOP/FEVAL

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isa(A, 'chebfun') )
    error('CHEBOP:mtimes:invalid', 'Operation is undefined.');
    
elseif ( isa(B, 'chebfun') )    % Let A operate on B
    if ( nargin(A.op) == 1 )
        C = feval(A, B);
    else
        C = feval(A, chebfun(@(x) x, A.domain), B);
    end
    
elseif ( isnumeric(A) || isnumeric(B) )
    % Switch argument to make sure A is numeric
    if isnumeric(B)
        temp = A; A = B; B = temp; %#ok<NASGU>
    end
    
    % Initialize C as a CHEBOP
    C = B;  % change this if ID's are put in chebops!
    
    % Multiplication of anonymous functions is not supported in Matlab. Need to
    % work around that.
    funString = func2str(C.op);                     % Anon. func. string
    firstRightParLoc = min(strfind(funString,')')); % First ) in string
    funArgs = funString(2:firstRightParLoc);        % Grab variables name
    C.op = eval(['@',funArgs,'A*C.op',funArgs]);    % Create new anon. func
    
elseif ( isa(A,'chebop') && isa(B,'chebop') )       
    % Chebop composition is not yet supported. It will probably be a mess once
    % we start getting systems involved, as the inputs need to be shuffled
    % correctly. A remedy might be through a nested function in this file.
    error('CHEBFUN:CHEBOP:MTIMES:COMPOSITIONS',...
        'Chebop composition is currently not supported.');   
end