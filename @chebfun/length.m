function [out, lengthFuns] = length(F)
%LENGTH   Length of a Chebfun.
%   LENGTH(F) returns the length of a scalar-valued CHEBFUN object F, which is
%   defined as the sum of the length of F.funs. If F is a quasimatrix, then
%   LENGTH(F) returns the maximum length of the columns.
%
%   [LEN, LENFUNS] = LENGTH(F) also returns the length of each of the piecewise
%   components of the scalar-valued CHEBFUN object F. If F is array-valued
%   LENFUNS = NaN.
%
% See also SIZE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) )
    out = 0;
else
    out = zeros(1, numel(F));
    for k = 1:numel(F)
        out(k) = sum(cellfun(@length, F(k).funs));
    end
    out = max(out);
end

if ( nargout > 1 && (numColumns(F) == 1) )
    lengthFuns = cellfun(@length, F(1).funs).';
    if ( F(1).isTransposed )
        lengthFuns = lengthFuns.';
    end
else
    lengthFuns = NaN;
end

end
