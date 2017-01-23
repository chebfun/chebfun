function out = isempty(f)
%ISEMPTY   True for an empty SINGFUN.
%   ISEMPTY(F) returns TRUE if F is an empty SINGFUN and FALSE otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
   
% Check if the smooth part is empty:
out = isempty(f.smoothPart);
    
end
