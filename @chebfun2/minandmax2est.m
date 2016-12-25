function varargout = minandmax2est(varargin)
%MINANDMAX2EST   Estimates the minimum and maximum of a CHEBFUN2.
%   mM = MINANDMAX2EST(F) returns estimates for the minimum and maximum of the
%   CHEBFUN2 F over its domain.  mM is a vector of length 2 such that mM(1) is
%   the estimated minimum and mM(2) is the estimated maximum.
%
%   mM = MINANDMAX2EST(F, N) returns estimates for the minimum and maximum of
%   the CHEBFUN2 F over its domain, based on samples on an N by N grid
%   (N = 33 by default).
% 
% See also MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = minandmax2est@separableApprox(varargin{:});

end
