function out = isempty( f ) 
%ISEMPTY   True for empty BALLFUN.
%   ISEMPTY(F) returns 1 if F is an empty BALLFUN and 0 otherwise. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isempty( f.coeffs ); 
end