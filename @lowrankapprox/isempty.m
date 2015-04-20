function out = isempty( f ) 
%ISEMPTY   True for empty LOWRANKAPPROX.
%   ISEMPTY(F) returns 1 if F is an empty LOWRANKAPPROX object and 0 otherwise. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isempty( f.cols ) && isempty( f.rows ); 

end
