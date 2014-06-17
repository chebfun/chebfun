function f = sum( f, dim )
%SUM   Definite Integration of a CHEBFUN2.
%   G = sum(F,DIM) where DIM is 1 or 2 integrates only over Y or X respectively,
%   and returns as its output a chebfun in the remaining variable.
%
%   G = sum(F) is the same as sum(F,1)
%
% See also SUM2. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    f = []; 
    return; 
end

% Default to y direction: 
if ( nargin == 1 )
    dim = 1;
end

% Get the low rank representation for f. 
[cols, D, rows] = cdr(f);
dom = f.domain; 

if ( dim == 1 )
    % Integrate over y: 
    f = rows * ( sum(cols) * D ).';
    if ( isa(f, 'chebfun') ) 
        f = simplify( f.' ); 
    else
        % f = double 
        f = chebfun(f, dom(1:2)).'; 
    end
elseif ( dim == 2 )
    f = cols * ( D * sum( rows ).' );
    if  ( isa(f, 'chebfun') ) 
        f = simplify( f );
    else
        % f = double 
        f = chebfun( f, dom(3:4) ); 
    end
else 
    error('CHEBFUN:CHEBFUN2:sum:unknown', ...
          'Undefined function ''sum'' for that dimension');
end

end
