function [x, w] = points(disc, kind)
%POINTS    Discretization points in COLLOC2.
%   X = POINTS(DISC) returns 2nd-kind points using the domain and dimension
%   stored in DISC.
%
%   [X, W] = POINTS(DISC) also returns Clenshaw-Curtis weights.
%
%   POINTS(DISC, KIND) overrides to allow returning 1st kind points if KIND=1. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Default to second kind points
if ( nargin < 2 )
    kind = 2;
end

% Obtain useful info
d = disc.domain;
numInt = disc.numIntervals;
n = disc.dimension;

% FIXME: This is deprecated usage. Is it necessary?  
if ( numel(n) == 1 )
    n = repmat(n, 1, numInt);
end

% Create output as cells for convenience. 
x = cell(numInt, 1);
w = cell(1, numInt);
for k = 1:numInt
    
    % To save time, don't call again unless the dimension has changed.
    if ( k==1 ) || ( n(k) ~= n(k-1) )
        if ( kind == 2 )
            [x0, w0] = chebtech2.chebpts(n(k));
        else
            [x0, w0] = chebtech1.chebpts(n(k));
        end
    end
    
    % The points and weights returned by the CHEBTECH methods above live on
    % [-1, 1]. Transform the points to the interval we're working on.
    dif = (d(k + 1) - d(k))/2;
    x{k} = x0*dif + (d(k +1) + d(k))/2;
    
    if ( nargout > 1 )
        w{k} = w0*dif;
    end
    
end

% Convert output to vectors. 
x = cell2mat(x);
if (nargout > 1)
    w = cell2mat(w);
end

end
