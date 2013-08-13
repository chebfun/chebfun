function f = transpose(f)
%.'   CHEBFUN transpose.
%   F.' is the non-conjugate transpose of F, i.e., it converts a column CHEBFUN
%   to a row CHEBFUN and vice versa.
%
% See also CHEBFUN/CTRANSPOSE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

f.isTransposed = ~f.isTransposed;

end
