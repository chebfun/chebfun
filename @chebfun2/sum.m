function f = sum( f, dim )
%SUM  Definite Integration of a chebfun2.
%
% G = sum(F,DIM) where DIM is 1 or 2 integrates only over Y or X
% respectively, and returns as its output a chebfun in the remaining
% variable.
%
% G = sum(F) is the same as sum(F,1)
%
% See also SUM2. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default to y direction: 
if ( nargin == 1 )
    dim = 1;
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues; 
d = 1./piv; 
d(d==inf) = 0;  % set infinite values to zero. 

if ( dim == 1 )
    % Integrate over y: 
    f = rows * ( sum(cols) * diag( d ) ).';
elseif ( dim == 2 )
    f = cols * ( diag( d ) * sum( rows ).' );
else 
    error('CHEBFUN2:SUM:unknown', ...
          'Undefined function ''sum'' for that dimension');
end
end