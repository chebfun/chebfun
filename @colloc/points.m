function [x, w, v] = points(varargin)
%POINTS    Discretization points.
%   X = COLLOC.POINTS(DISC, KIND) returns KIND-kind points using the domain and
%   dimension stored in DISC. KIND must be either 1 or 2.
%
%   An alternate calling sequence is COLLOC.POINTS(DOMAIN, DIMENSION, KIND).
%
%   [X, W, V] = COLLOC.POINTS(DISC, KIND) also returns Clenshaw-Curtis
%   weights and barycentric interpolation weights, respectively.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: I find it weird that this exists. NH Apr 2014

if ( nargin == 2 )
    disc = varargin{1};
    kind = varargin{2};
    d = disc.domain;
    n = disc.dimension;
    
elseif ( nargin == 3 )
    d = varargin{1};
    n = varargin{2};
    kind = varargin{3};
    
else
    error('Must be called with two or three arguments.')
    
end

numInt = length(n);

% Create output as cells for convenience.
x = cell(numInt, 1);
w = cell(1, numInt);
v = cell(numInt, 1);
for k = 1:numInt
    
    % To save time, don't call again unless the dimension has changed.
    if ( k==1 ) || ( n(k) ~= n(k-1) )
        if ( kind == 2 )
            [x0, w0, v0] = chebtech2.chebpts(n(k));
        else
            [x0, w0, v0] = chebtech1.chebpts(n(k));
        end
    end
    
    % The points and weights returned by the CHEBTECH methods above live on
    % [-1, 1]. Transform the points to the interval we're working on.
    dif = (d(k + 1) - d(k))/2;
    x{k} = x0*dif + (d(k +1) + d(k))/2;
    
    if ( nargout > 1 )
        w{k} = w0*dif;
        if ( nargout > 2 )
            v{k} = v0;
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
end

end