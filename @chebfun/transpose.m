function f = transpose(f)
%.'   Transpose.
%   F.' is the non-conjugate transpose of F, i.e., it converts a column chebfun
%   to a row chebfun and vice versa.
%
% See also CTRANSPOSE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

f.isTransposed = ~f.isTransposed;

end
