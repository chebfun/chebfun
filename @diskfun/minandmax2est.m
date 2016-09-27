function varargout = minandmax2est(varargin)
%MINANDMAX2EST   Estimates the minimum and maximum of a DISKFUN.
%   mM = MINANDMAX2EST(F) returns estimates for the minimum and maximum of the
%   DISKFUN F over its domain.  mM is a vector of length 2 such that mM(1) is
%   the estimated minimum and mM(2) is the estimated maximum.
%
%   mM = MINANDMAX2EST(F, N) returns estimates for the minimum and maximum of
%   the DISKFUN F over its domain, based on the evaluation on an N by N
%   Chebyshev grid in the domain of F (N = 33 by default).
%
% See also DISKFUN/MINANDMAX2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f = varargin{1};
if ( isempty(f) )
    mM = [0, 0];
    varargout = {mM};
    return
end

% Give CDR to separableApprox so it evaluates in polar.
f = cart2pol(f, 'cdr');
varargin{1} = f;

varargout{1:nargout} = minandmax2est@separableApprox(varargin{1:nargin});

end
