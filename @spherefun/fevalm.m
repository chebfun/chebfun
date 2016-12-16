function varargout = fevalm(varargin)
% FEVALM   Evaluate a SPHEREFUN in spherical coordinates.
% 
% Z = FEVALM(F, LAM, TH) returns a matrix of values Z of size
% length(TH)-by-length(LAM). (LAM,TH) are the spherical coordinates for the
% evaluation points on the sphere, with -pi <= LAM <= pi the azimuthal
% angle and 0 <= TH <= pi the elevation (polar) angle (both measured in
% radians). LAM and TH should be vectors of doubles. This is equivalent to
% making a meshgrid of the vectors LAM and TH and then using FEVAL to
% evaluate at that grid.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = fevalm@separableApprox(varargin{:});

end
