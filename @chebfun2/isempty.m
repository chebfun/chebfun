function bol = isempty( f ) 
% ISEMPTY True for empty chebfun2
%
% ISEMPTY(F) returns 1 if F is an empty chebfun2 object and 0 otherwise. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

bol_cols = isempty( f.cols ); 
bol_rows = isempty( f.rows ); 
bol = bol_cols & bol_rows; 

end