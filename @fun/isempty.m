function out = isempty(f)
%ISEMPTY   True for an empty FUN.
%   ISEMPTY(F) returns TRUE if F is an empty FUN and FALSE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the ONEFUN is empty:
out = isempty(f.onefun);


end
