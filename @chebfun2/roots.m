function varargout = roots(varargin)
%ROOTS   Zero contours of a CHEBFUN2.
%   R = ROOTS(F), returns the zero contours of F as a quasimatrix of chebfuns.
%   Each column of R is one zero contour. This command only finds contours when
%   there is a change of sign and it can also group intersecting contours in a
%   non-optimal way. Contours are computed to, roughly, four digits of
%   precision. In particular, this command cannot reliably compute isolated real
%   roots of F or zero curves lying close to the boundary of the domain.
%
%   In the special case when F is of length 1 then the zero contours are found
%   to full precision.
%
%   R = ROOTS(F, G) returns the isolated points of F and G.
%
%   R = ROOTS(F, G, METHOD) the underlying rootfinding algorithm can be
%   supplied. If METHOD = 'ms' or METHOD = 'marchingsquares', then the Marching
%   Squares algorithm is employed. The Marching Squares algorithm is fast but
%   not particularly robust. If METHOD = 'resultant', then a hidden variable
%   resultant method based on Bezout resultants is employed. The Resultant
%   method is slower but far more robust. See the CHEBFUN2V/ROOTS()
%   documentation to see which algorithm is used when no METHOD is passed.
%
% See also CHEBFUN2V/ROOTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = roots@separableApprox(varargin{:});

end
