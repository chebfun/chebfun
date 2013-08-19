function out = length(f)
%LENGTH   Length of a Chebfun.
%  LENGTH(F) returns the length of a CHEBFUN object F, which is defined as the
%  sum of the length of f.funs.
%
% See also CHEBFUN/SIZE.

out = sum(cellfun(@length, f.funs));

end