function vals = feval(f, varargin)
%FEVAL   Evaluate a BALLFUNV at one or more points.
%   Y = FEVAL( F, X, Y, Z) evaluates a BALLFUNV F at a point (X,Y,Z) in Cartesian
%   coordinates, where X, Y and Z are doubles.
%
%   Y = FEVAL( F, R, LAM, TH, 'spherical') evaluates a ballfun F in
%   spherical coordinates (R,LAM,TH). Here R, LAM and THETA are doubles representing 
%   the radius, azimuthal and polar angles (in radians) and must be points 
%   in the unit ball.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
    vals = [];
    return
end

% Extract the components
F = f.comp;

% Evaluate the ballfun objects
V1 = feval(F{1}, varargin{:});
V2 = feval(F{2}, varargin{:});
V3 = feval(F{3}, varargin{:});

% Return the values
vals = [V1, V2, V3];
end
