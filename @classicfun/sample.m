function [points, values] = sample(f, varargin)
%SAMPLE   Sample a CLASSICFUN on an "appropriate" grid.
%   [POINTS, VALUES] = SAMPLE(F, N) returns a vector POINTS of N
%   "appropriately-chosen" points and a vector VALUES of the values of F at
%   those points.  What "appropriately-chosen" means depends on the type of
%   representation on which F is ultimately based.
%
%   [POINTS, VALUES] = SAMPLE(F) uses N = LENGTH(F).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [points, values] = sample(f.onefun, varargin{:});
    point = f.mapping.For(points);
end
