function out = isempty(f)
%ISEMPTY   True for an empty SINGFUN.
%   ISEMPTY(F) returns TRUE if F is an empty SINGFUN and FALSE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    % An array cannot be empty.
    out = false;
elseif ( numel(f) == 1 )
    % Check if the smooth part is empty:
    out = isempty(f.smoothPart);
else 
    % numel(f) == 0, so f must be empty.
    out = true;
end

end
