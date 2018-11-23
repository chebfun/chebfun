function b = iszero(f)
%ISZERO   Check if a BALLFUNV is identically zero on its domain.
%   OUT = ISZERO( F ) return 1 if the BALLFUNV is the zero vector up to 
%   machine precision, and 0 otherwise.
%
%   See also ISEQUAL.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
b = (iszero(F{1}) && iszero(F{2}) && iszero(F{3}));

end
