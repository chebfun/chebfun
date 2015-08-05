function G = divgrad(F)
%DIVGRAD   Laplacian of a CHEBFUN2V.
%   F = DIVGRAD(F) returns the Laplacian of a CHEBFUN2V i.e.,
%       divgrad(F) = F(1)_xx + F(2)_yy
%
% This command is not defined for a chebfun2v with 3 components. 
%
% Also see CHEBFUN2V/LAP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

nComponents = F.nComponents; 
if ( nComponents > 2 ) 
    error('CHEBFUN:CHEBFUN2V:divgrad:components',...
        'Command is not defined for CHEBFUN2V objects with >2 components.')
end
     
Fc = F.components; 
G = diff(Fc{1}, [2,0]) + diff(Fc{2}, [0,2]);
 
end
