function H = mldivide(F, G)
%\   Left divide for CHEBFUN3V objects.
%   c\F Divides each component of a CHEBFUN3V by the scalar c.
%
%   Only allowed to divide a CHEBFUN3V by a scalar.
%
% See also CHEBFUN3V/MRDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(G) ) )
    H = chebfun3v;
    return
end

if ( ~isa(F, 'double') )
    error('CHEBFUN:CHEBFUN2V:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide.
H = G;
for j = 1:G.nComponents
    H.components{j} = mldivide(F, G.components{j});
end

end