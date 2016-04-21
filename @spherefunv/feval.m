function vals = feval(F, varargin)
%FEVAL   Evaluate a SPHEREFUNV at one or more points.
%   Y = FEVAL(F, LAMBDA, THETA) evaluates F at (LAMBDA, THETA) where
%   LAMBDA and THETA are doubles representing the longitudinal (or
%   azimuthal) and latitudinal (or elevation) angles.
%
%   Y = FEVAL(F, X, Y, Z)  evaluates F at a point (X,Y,Z) in
%   Cartesian cooridnates on the surface of a sphere.  
%
% See also SPHEREFUNV/SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) ) 
    vals = []; 
    return
end

if nargin < 3
    error('SPHEREFUN:SPHEREFUNV:feval:numInputs', ...
        'Wrong number of evaluation points.');
end

vals = zeros(3, numel(varargin{1})); 

% Evaluate each component:
for jj = 1:3
   vals(jj, :) = feval(F.components{jj}, varargin{:});  
end

end
