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

if ( ~isa(F, 'double') && ~isa(F, 'chebfun3') )
    error('CHEBFUN:CHEBFUN3V:mldivide:nonScalar', ...
        'Division must be by a scalar or a CHEBFUN3.');
end

H = G; 
% Left Componentwise divide:
if ( isa(F,'double') )
    for j =1:G.nComponents
        H.components{j} = mldivide(F, G.components{j});
    end
else
    for j = 1:G.nComponents
        H.components{j} = ldivide(F, G.components{j});
    end
end
    
end