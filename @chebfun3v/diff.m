function F = diff(F, k, dim)
%DIFF   Componentwise partial derivative of a CHEBFUN3V object.
%   DIFF(F) is the derivative of components of F along the first variable.
%
%   DIFF(F, K) is the Kth derivative of each component of F along the first
%   variabel.
%
%   DIFF(F, K, DIM) is the Kth derivative of F along the dimension DIM.
%   DIM = 1 (default) is the derivative in the 1st input variable.
%   DIM = 2 is the derivative in the 2nd variable.
%   DIM = 3 is the derivative in the 3rd varialbe.
%
%   DIFF(F, [K1 K2], [DIM1 DIM2]) means K1-th derivative of F in dimension
%   DIM1 and K2-th derivative in dimension DIM2. DIM1 and DIM2 can be 1, 2 
%   or 3 in any order.
%
%   DIFF(F, [K1 K2 K3]) is the K1-th partial derivative of F in the first 
%   variable, K2-th partial derivative of F in the second variable and 
%   K3-th partial derivative of F in the third variable. 
%   For example, DIFF(F, [1 2 3]) is d^6F/(dx d^2y d^3z).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) ) 
    return
end

% Defaults:
if ( ( nargin == 1 ) || isempty(k) )
    k = 1;
end
if ( nargin < 3 ) 
    dim = 1; 
end

% Diff each component. 
for j = 1:F.nComponents
    F.components{j} = diff(F.components{j}, k, dim); 
end

end