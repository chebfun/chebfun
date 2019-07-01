function b = iszero( f )
%ISZERO   Check if a BALLFUN is identically zero on its domain.
%   ISZERO( F ) returns 1 if the BALLFUN is exactly the zero function,
%   and 0 otherwise. 
%
%   See also ISEQUAL.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
b = nnz(F) == 0;
end
