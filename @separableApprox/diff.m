function F = diff(F, k, dim)
%DIFF   Derivative of a SEPARABLEAPPROX object.
%   DIFF(F) is the derivative of F along the y direction.
%
%   DIFF(F, N) is the Nth derivative of F in the y direction.
%
%   DIFF(F, N, DIM) is the Nth derivative of F along the dimension DIM.
%     DIM = 1 (default) is the derivative in the y-direction.
%     DIM = 2 is the derivative in the x-direction.
%
%   DIFF(F, [NX NY]) is the partial derivative of NX of F in the first variable,
%   and NY of F in the second derivative. For example, DIFF(F,[1 2]) is
%   d^3F/dxd^2y.
%
% See also GRADIENT, SUM, PROD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    return
end

% Default to first derivative:
if ( nargin < 2 )
    k = 1;
end

% Default to partial derivative in y:
if ( nargin < 3 )
    dim = 1;
elseif ( numel(dim) ~= 1 )
    error('CHEBFUN:SEPARABLEAPPROX:diff:dim1', 'Dim should be either 1 or 2.');
end

% Diff the individual column and row slices:
if ( numel(k) == 2 && nargin < 3 )
   F.cols =  diff(F.cols, k(2));
   F.rows = diff(F.rows, k(1));
elseif ( dim == 1 )
    F.cols = diff(F.cols, k);
elseif ( dim == 2 )
    F.rows = diff(F.rows, k);
else 
    error('CHEBFUN:SEPARABLEAPPROX:diff:dim2', ...
        'Can compute derivative in x or y only.');
end

end
