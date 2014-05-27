function f = diff(f, varargin) %diff(f,order,dim,varargin)
%DIFF Derivative of a ADchebfun2.
%
% DIFF(F) is the derivative of F along the y direction.
%
% DIFF(F,N) is the Nth derivative of F in the y direction.
%
% DIFF(F,N,DIM) is the Nth derivative of F along the dimension DIM.
%     DIM = 1 (default) is the derivative in the y-direction.
%     DIM = 2 is the derivative in the x-direction.
%
% DIFF(F,[NX NY]) is the partial derivative of NX of F in the first
% variable, and NY of F in the second derivative. For example, DIFF(F,[1
% 2]) is d^3F/dxd^2y.
%
% See also GRADIENT, SUM, PROD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


% Start by differentiating the chebfun2 of the ADchebfun
f.chebfun2 = diff(f.chebfun2,varargin{:});

% Find out how many derivatives are required in each direction

if ( nargin == 1 ) % defaults.
    nx = 0;
    ny = 1;
elseif ( nargin == 2 )
    % Two arguments passed, so second argument must be the order
    order = varargin{1};
    if length(order) == 1 % diff in y is default.
        nx = 0;
        ny = order;
    elseif length(order) == 2
        % Got a vector passed, first entry is number of derivatives in x
        % direction, second entry is number of derivatives in y direction.
        nx = order(1);
        ny = order(2);
    else
        error('CHEBFUN2:DIFF','Undetermined direction of differentiation.');
    end
else
    % Three arguments passed, second one is the order, third is the dimension
    order = varargin{1};
    dim = varargin{2};
    
    if dim == 1
        nx = 0;
        ny = order;
    else
        nx = order;
        ny = 0;
    end    
end

% Update the derivative information
f.der = diff(f.der, nx, ny);

end