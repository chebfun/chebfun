function varargout = fevalm(varargin)
%FEVALM   Evaluate a DISKFUN in polar coordinates.
% 
%   Z = FEVALM(F, THETA, R) returns a matrix of values Z of size
%   length(R)-by-length(THETA). (R,THETA) are polar coordinates for the
%   evaluation points on the disk.  They should be vectors of doubles. 
%   Calling this function as above is equivalent to making a meshgrid of
%   the vectors THETA and R and then using FEVAL to evaluate at that grid.
%
% See also DISKFUN/FEVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = fevalm@separableApprox(varargin{:});

end
