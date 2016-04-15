function f = divergence( F )
%DIVERGENCE   Divergence of a CHEBFUN3V.
%   DIVERGENCE(F) returns the divergence of the CHEBFUN3V i.e.,
%       divergence(F) = F_x + F_y + F_z.
%
% See also CHEBFUN3V/DIV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    f = chebfun3;
    return
end

Fc = F.components; 
diff1 = diff(Fc{1}, 1, 1);
diff2 = diff(Fc{2}, 1, 2);
diff3 = diff(Fc{3}, 1, 3);

f = diff1 + diff2 + diff3; % Use compressed plus.

end