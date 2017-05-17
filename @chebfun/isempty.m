function out = isempty(f)
%ISEMPTY   Test for empty CHEBFUN.
%   ISEMPTY(F) returns logical true if F is an empty CHEBFUN and false
%   otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If there is a .fun and f.funs(1) is not empty the f is not empty.
if ( numel(f) > 1 || ( (numel(f.funs) > 0) && ~isempty(f.funs{1}) ) )
    out = false;
else
    out = true;
end

end
