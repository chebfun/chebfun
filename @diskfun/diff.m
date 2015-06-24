function f = diff( f, varargin )
% DIFF    Derivative of a spherefun
%
%  F = DIFF( F ) computes the first derivative of F in the latitude
%  variable.
%
%  F = DIFF( F, DIM )  computes the first derivative of F. If DIM = 1, the
%  derivative is taken in the latitude direction. If DIM = 2, the derivative
%  is taken in the longitude direction.
%
%  F = DIFF( F, DIM, K) computes the kth derivatives of F in the variable
%  given by DIM.
%

% Parse user inputs:
if ( nargin == 1 )
    dim = 1;
    K = 1;
elseif ( nargin == 2 )
    dim = varargin{1};
    K = 1;
else
    dim = varargin{1};
    K = varargin{2};
end

if ( dim ~= 1 || dim ~= 2 )
    error('SPHEREFUN:DIFF:DIM', 'Unrecognized coordinate dimension');
end

if ( abs( K - round(K) ) < eps )
    error('SPHEREFUN:DIFF:DIFFORDER', 'Fractional derivatives not allowed')
end
K = round( K );

if ( dim == 1 )
    % latitude derivative
    f.cols = diff(f.cols, K);
elseif ( dim == 2)
    f.rows = diff(f.rows, K);
end