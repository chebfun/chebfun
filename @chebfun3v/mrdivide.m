function H = mrdivide(F, G)
%/   Right divide for CHEBFUN3V objects.
%   F/c divides each component of a CHEBFUN3V by a scalar.
%
%   Only allowed to divide by scalars.
% 
% See also CHEBFUN3V/MLDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) || isempty(G) )
   H = chebfun3v();
   return 
end

if ( ~isa(G, 'double') && ~isa(G, 'chebfun3') )
    error('CHEBFUN:CHEBFUN3V:mrdivide:nonScalar', ...
        'Division must be by a scalar or a CHEBFUN3 object.');
end

H = F; 
% Componentwise divide:
if ( isa(G,'double') )
    for j =1:F.nComponents
        H.components{j} = mrdivide(F.components{j}, G);
    end
else
    for j = 1:F.nComponents
        H.components{j} = rdivide(F.components{j}, G);
    end
end

end