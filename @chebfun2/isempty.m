function out = isempty( f ) 
%ISEMPTY   True for empty CHEBFUN2.
%   ISEMPTY(F) returns 1 if F is an empty CHEBFUN2 object and 0 otherwise. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

out = isempty( f.cols ) && isempty( f.rows ); 

end