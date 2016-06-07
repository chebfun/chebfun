function out = isempty( f ) 
%ISEMPTY   True for empty SEPARABLEAPPROX.
%   ISEMPTY(F) returns 1 if F is an empty SEPARABLEAPPROX object and 0 otherwise. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isempty( f.cols ) && isempty( f.rows ); 

end
