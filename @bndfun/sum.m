function out = sum(f, dim)
%SUM   Definite integral of a BNDFUN on its interval [a,b]
%   S = SUM(F) is the integral from a to b of F.
%
%   If F is an array-valued BNDFUN then the result is a row vector containing
%   the definite integrals of each column.
%
%   SUM(F, 2) sums over the second dimension of F, i.e., adds up its columns. If
%   F is a scalar-valued BNDFUN, this simply returns F.
%
% See also CUMSUM, DIFF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Sum across array-valued BNDFUN columns if dim = 2:
if ( (nargin > 1) && (dim == 2) )
    f.onefun = sum(f.onefun, dim);
    out = f;
    return
end

%%
% Otherwise, compute the sum column-wise.

% Rescaling factor, (b - a)/2
rescaleFactor = 0.5*diff(f.domain);

% Assign the output to be the sum of the onefun of the input, rescaled.
out = sum(f.onefun)*rescaleFactor;

end
