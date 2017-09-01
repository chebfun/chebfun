function varargout = median(varargin)
%MEDIAN       Median value of a CHEBFUN2.
%   G = MEDIAN(F) returns a CHEBFUN G representing the median of the CHEBFUN2
%   along the y direction, i.e., G = @(x) median( F ( x, : ) ).
%
%   G = MEDIAN(F, DIM) returns a CHEBFUN G representing the median of F along
%   the direction given by DIM, i.e., y-direction if DIM = 1 and x-direction if
%   DIM = 2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = median@separableApprox(varargin{:});

end
