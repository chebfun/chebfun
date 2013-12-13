function f = imag(f)
%IMAG   Complex imaginary part of a CHEBFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Handle the empty case:
if ( isempty(f) )
    return
end

% Take imaginary part of the pointValues:
f.pointValues = imag(f.pointValues);

% Take imaginary part of the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = imag(f.funs{k});
end

end
