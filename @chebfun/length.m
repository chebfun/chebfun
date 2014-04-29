function out = length(F)
%LENGTH   Length of a Chebfun.
%   LENGTH(F) returns the length of a scalar-valued CHEBFUN object F, which is
%   defined as the sum of the length of F.funs. If F is an quasimatrix, then
%   LENGTH(F) returns the maximum length of the columns.
%
% See also SIZE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
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

end
