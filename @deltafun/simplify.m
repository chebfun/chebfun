function f = simplify(f, varargin)
%SIMPLIFY  Simplifys a DELTAFUN object F.
%
% See also SUM, CUMSUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%TODO: NH: DOCS

f.funPart = simplify(f.funPart, varargin{:});
f = simplifyDeltas(f, varargin{:});

end