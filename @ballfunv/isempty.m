function out = isempty(F)
%ISEMPTY   Empty boolean check for a BALLFUNV object. 
%   ISEMPTY(F) returns 1 if every component of F is an empty BALLFUNV, and
%   returns 0 otherwise.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( F.comp) ) 
    out = 1; 
    return
end

% Take isempty of each component:
out = cellfun(@isempty, F.comp, 'UniformOutput', false);
out = all(cell2mat(out));

end