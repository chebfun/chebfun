function b = iszero(f)
%ISZERO   Check if a BALLFUNV is identically zero on its domain.
%   ISZERO( F ) returns 1 if the BALLFUNV is the zero vector and 0 
%   otherwise.
%
%   See also ISEQUAL.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
    b = 1;
    return
end

F = f.comp;
b = (iszero(F{1}) && iszero(F{2}) && iszero(F{3}));
end
