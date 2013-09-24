function out = innerProduct(f, g)
%INNERPRODUCT   Compute the inner product of two CHEBFUN objects.
%   INNERPRODUCT(F, G) returns the L2 inner product of the two CHEBFUN objects F
%   and G (conjugate linear in F).
%
%   If F and/or G are array-valued CHEBFUN objects, then the result is a matrix
%   whose i,j entry is the inner product of the ith column of F with the jth
%   column of G.
%
% See also NORM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Overlap the CHEBFUN objects:
[f, g] = overlap(f, g);

% Initialise the output:
out = zeros(min(size(f)), min(size(g)));

% Loop over the FUNs:
for k = 1:numel(f.funs)
    out = out + innerProduct(f.funs{k}, g.funs{k});
end 

end
