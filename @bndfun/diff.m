function f = diff(f,k,varargin)
%IMAG	Derivative of a BNDFUN.
%   DIFF(F) is the derivative of F and DIFF(F, K) is the K-th derivative.
%
%   DIFF(F, K, DIM), where DIM is one of 1 or 2, takes the Kth difference
%   along dimension DIM. For DIM = 1, this is the same as above. For DIM =
%   2, this is a finite difference along the columns of a vectorised
%   FUNCHEB. If F has L columns, a BNDFUN with an empty onefun property
%   will be returned for K >= L.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

if nargin < 2
    k = 1;
end

% Rescaling factor, (b-a)/2, to the kth power
dab05k = (.5*diff(f.domain))^k;

% Assign the onefun of the output to be the derivative of the onefun of the
% input, and scale the result:
f.onefun = diff(f.onefun,k,varargin{:})/dab05k;

end
