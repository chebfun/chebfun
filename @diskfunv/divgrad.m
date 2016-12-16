function G = divgrad(F)
%DIVGRAD   Laplacian of a DISKFUNV.
%   F = DIVGRAD(F) returns the Laplacian of a DISKFUNV i.e.,
%       divgrad(F) = F(1)_xx + F(2)_yy
%
% Also see DISKFUNV/LAPLACIAN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) )
    G = diskfun();
    return
end

% Extract components: 
Fc = F.components; 

% Calculate DIVGRAD:
G = diff(Fc{1},1,2) + diff(Fc{2},2,2);
 
end