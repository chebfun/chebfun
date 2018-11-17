function b = iszero(f)
% ISZERO   Check if a BALLFUNV is identically zero on its domain.
%   ISZERO(F) returns true if F is identically zero or empty on F.domain and
%   false otherwise.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
b = (iszero(F{1}) && iszero(F{2}) && iszero(F{3}));

end
