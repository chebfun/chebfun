function varargout = median(varargin)
%MEDIAN   Median value of a DISKFUN.
%   G = MEDIAN(F) returns a CHEBFUN G representing the median of the DISKFUN
%   along the radial direction, i.e., G = @(t) median( F ( t, : ) ), where
%   F(theta, r) expresses F in polar coordinates and -pi <= theta <= pi, 
%   0 <= r <= 1. 
%
%   G = MEDIAN(F, DIM) returns a CHEBFUN G representing the median of F along
%   the direction given by DIM, i.e., radial direction if DIM = 1 and angular 
%   direction if DIM = 2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = varargin{1}; 
if isempty(f)
    varargout ={[]};
    return
end
f = cart2pol(f); %gives back polar f restricted to [-pi pi 0 1]
[varargout{1:nargout}] = median@separableApprox(f,varargin{2:end});

end