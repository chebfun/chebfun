function varargout = roots(varargin)
%ROOTS   Zero contours of a CHEBFUN2.
%   R = ROOTS(F) returns the zero contours of F as a quasimatrix of chebfuns.
%   Each column of R is one zero contour. This command only finds contours when
%   there is a change of sign and it can also group intersecting contours in a
%   non-optimal way.
%
%   For a faster plot to graphical accuracy use CONTOUR(F, [0 0]).
%
%   R = ROOTS(F, G, METHOD) allows the underlying rootfinding algorithm to
%   be specified.  If METHOD = 'ms' or 'marchingsquares', the Marching
%   Squares algorithm is employed, which is fast but not very robust.
%   If METHOD = 'resultant', a hidden variable resultant method
%   based on Bezout resultants is employed, slower but more robust.
%   See the CHEBFUN2V/ROOTS documentation to see which algorithm is used
%   when no METHOD is passed.
%
% Example:
%   cheb.xy;
%   f = x.^2 + y.^2 - 1/4;
%   roots(x,f)                               % [0 -.5; 0 .5]
%   c = roots(f);
%   arclength = sum(abs(diff(c)))            % pi
%   area = abs(sum(real(c).*diff(imag(c))))  % pi/4
%
% See also CHEBFUN2V/ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = roots@separableApprox(varargin{:});

end
