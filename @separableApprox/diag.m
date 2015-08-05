function f = diag( f, varargin )
%DIAG(F)   Diagonal of a SEPARABLEAPPROX.
%   G = DIAG(F) returns the CHEBFUN representing g(x) = f(x, x).
%
%   G = diag(F,C) returns the CHEBFUN representing g(x) = f(x, x+c).
%
% See also TRACE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty( f ) ) 
    f = chebfun;
    return
end 

if ( nargin == 1 )     
    % Default to zero diagonal shift.
    c = 0; 
else 
    c = varargin{1};
    if ( ~isa(c, 'double') )
        error('CHEBFUN:SEPARABLEAPPROX:diag:diag', ...
            'Second argument to diag should be a double.');
    end
end

dom = f.domain;
dom = [max(dom(1), dom(3)-c), min(dom(2), dom(4)-c)]; % Find domain of diagonal.
f = chebfun( @(x) feval( f, x, x + c ), dom );        % Construct the diagonal. 

end
