function f = simplify(f, varargin)
%SIMPLIFY  Simplifys a DELTAFUN object F.
%   F = SIMPLIFY(F) simplifies the FUNPART and removes any trivial delta
%   functions from the delta fun F. Preferences can be passed via VARARGIN.
%
% See also SIMPLIFY, CLEANROWS, CLEANCOLUMNS

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

f.funPart = simplify(f.funPart, varargin{:});
f = simplifyDeltas(f, varargin{:});

end
