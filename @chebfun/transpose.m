function F = transpose(F)
%.'   CHEBFUN transpose.
%   F.' is the non-conjugate transpose of F, i.e., it converts a column CHEBFUN
%   to a row CHEBFUN and vice versa.
%
% See also CTRANSPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for k = 1:numel(F)
    F(k).isTransposed = ~F(k).isTransposed;
end

end
