function G = divgrad(F)
%DIVGRAD   Laplacian of a CHEBFUN3V object.
%   F = DIVGRAD(F) returns the Laplacian of a CHEBFUN3V i.e.,
%                                divgrad(F) = F(1)_xx + F(2)_yy + F(3)_zz.
%
%   This command is defined only for a CHEBFUN3V with 3 components.
%
% Also see CHEBFUN3V/LAP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

nComponents = F.nComponents; 
if ( (nComponents < 3) || (nComponents > 3) ) 
    error('CHEBFUN:CHEBFUN3V:divgrad:components',...
        'Defined only for CHEBFUN3V objects with three components.')
end
     
Fc = F.components; 
G = diff(Fc{1}, 2, 1) + diff(Fc{2}, 2, 2) + diff(Fc{3}, 2, 3); 
 
end