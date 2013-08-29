function out = isempty(f)
%ISEMPTY   True for an empty CHEBTECH.
%   ISEMPTY(F) returns TRUE if F is an empty CHEBTECH and FALSE otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( numel(f) > 1 )
    % An array cannot be empty.
    out = false;
elseif ( numel(f) == 1 )
    % Check if the values are empty:
    out = isempty(f.values);
else 
    % numel(f) == 0, so f must be empty.
    out = true;
end

end
