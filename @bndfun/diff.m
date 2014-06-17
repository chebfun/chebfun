function f = diff(f, k, dim)
%DIFF   Derivative of a BNDFUN.
%   DIFF(F) is the derivative of F, and DIFF(F, K) is the Kth derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference along
%   dimension DIM. For DIM = 1, this is the same as above. For DIM = 2, this is
%   a finite difference along the columns of an array-valued BNDFUN. If F has L
%   columns, a BNDFUN with an empty ONEFUN property will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( (nargin < 2) || isempty(k) )
    % Assume first derivative by default
    k = 1;
end

if ( nargin < 3 )
    % Assume dim = 1 by default
    dim = 1;
end

% Rescaling factor, (b-a)/2, to the kth power
rescaleFactork = (.5*diff(f.domain))^k;

% Assign the ONEFUN of the output to be the output of the DIFF method of the
% ONEFUN of the input. If we called DIFF with third argument equal to 2 (i.e.
% dim = 2), we only wanted to compute difference between columns, in which case,
% we should not rescale the result.
if ( dim == 1 )
    f.onefun = diff(f.onefun, k, dim)/rescaleFactork;
else
    f.onefun = diff(f.onefun, k, dim);
end

end
