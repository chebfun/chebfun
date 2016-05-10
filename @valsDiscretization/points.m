function [x, w, v, t] = points(varargin)
%POINTS    Discretization points.
%   X = VALSDISCRETIZATION.POINTS(DISC, POINTSFUN) returns points to be for the 
%   domain and dimension stored in DISC. POINTSFUN should be a function handle 
%   taking a single scalar input argument N which returns the required points in 
%   the interval [-1,1].
%
%   VALSDISCRETIZATION.POINTS(DOMAIN, DIMENSION, POINTSFUN) is an equivalent 
%   syntax.
%
%   [X, W, V] = VALSDISCRETIZATION.POINTS(DISC, POINTSFUN) also returns 
%   quadrature weights and barycentric interpolation weights, respectively, 
%   appropriately scaled to DISC.DOMAIN, if these are also returned by 
%   POINTSFUN.
%
%   [X, W, V, T] = VALSDISCRETIZATION.POINTS(DISC, POINTSFUN) returns also the 
%   angles T so that T = COS(X). This can often be computed more accurately than
%   X itself.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 )
    disc = varargin{1};
    pointsFun = varargin{2};
    d = disc.domain;
    n = disc.dimension;
    
elseif ( nargin == 3 )
    d = varargin{1};
    n = varargin{2};
    pointsFun = varargin{3};
    
else
    error('CHEBFUN:VALSDISCRETIZATION:points:nargin', ...
        'Must be called with two or three arguments.')
    
end

numInt = length(n);

% Create output as cells for convenience.
x = cell(numInt, 1);
w = cell(1, numInt);
v = cell(numInt, 1);
t = cell(numInt, 1);

for k = 1:numInt
    
    % To save time, don't call again unless the dimension has changed.
    if ( k == 1 ) || ( n(k) ~= n(k-1) )
        [varargout{1:nargout}] = pointsFun(n(k));
        x0 = varargout{1};
    end
    
    % The points and weights returned by the CHEBTECH methods above live on
    % [-1, 1]. Transform the points to the interval we're working on.
    x{k} = d(k+1)*(x0 + 1)/2 + d(k)*(1 - x0)/2;
    
    if ( nargout > 1 )
        dif = (d(k+1) - d(k))/2;
        w{k} = varargout{2}*dif;
        if ( nargout > 2 )
            v{k} = varargout{3};
            if ( nargout > 3 )
                t{k} = varargout{4};
            end
        end
    end
    
end

% Convert output to vectors.
x = cell2mat(x);
if ( nargout > 1 )
    w = cell2mat(w);
    if ( nargout > 2 )
        v = cell2mat(v);
    end
    if ( nargout > 3 )
        t = cell2mat(t);
    end
end

end
