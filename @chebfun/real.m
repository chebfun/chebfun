function f = real(f)
%REAL   Complex real part of a CHEBFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Handle the empty case:
if ( isempty(f) )
    return
end

% Take real part of the pointValues:
f.pointValues = real(f.pointValues);

% Take real part of the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = real(f.funs{k});
end

end
