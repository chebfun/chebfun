function f = sign(f, varargin)
%SIGN   Signum of a SIGN object.
%   SIGN(F) returns the sign of F, where F is a FUN object with no roots in
%   F.domain. If ~isempty(roots(F)), then SIGN(F) will return garbage with no
%   warning.
%
%   For the nonzero elements of complex F, sign(F) = F ./ ABS(F).
%
% See also ABS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Update preferences:
if ( nargin > 1 )
    pref = varargin{1};
    varargin{1} = f.onefun.pref(pref, pref.fun);
end

% Take the absolute value of the ONEFUN:
f.onefun = sign(f.onefun, varargin{:});

end