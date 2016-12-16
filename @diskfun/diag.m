function F = diag(f,varargin)
%DIAG(F)   Diagonal of a DISKFUN.
%   G = DIAG(F) returns the CHEBFUN representing g(r) = f(pi/4, r), where r
%   is the radial variable and is taken across the whole disk: -1 <= r <= 1
%
%   G = diag(F,T) returns the CHEBFUN representing g(r) = f(T, r), where
%  -pi < T < pi and r is the radial variable.
%
% See also DISKFUN/TRACE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )
    F = chebfun;
    return
end

% When c = 0, choose the diagonal radial slice for t = pi/4
if (nargin<2)
    F = feval(f, pi/4,':', 'polar');
else
    c = varargin{1};
    if abs(c) >= pi
        error('CHEBFUN:DISKFUN:diag:angleOutOfRange','Diagonal parameter must be between -pi and pi.');
    else
        F = feval(f, varargin{1},':', 'polar');
    end
end

end