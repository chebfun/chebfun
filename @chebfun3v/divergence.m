function f = divergence(F)
%DIVERGENCE   Divergence of a CHEBFUN3V object.
%   DIVERGENCE(F) returns divergence of the CHEBFUN3V object F as a 
%   CHEBFUN3. If F = U i + V j + W k, then divergence(F) = U_x + V_y + W_z.
%
% See also CHEBFUN3V/DIV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    f = chebfun3();
    return
end


if ( F.nComponents == 3 )
    Fc = F.components; 
    diff1 = diff(Fc{1}, 1, 1);
    diff2 = diff(Fc{2}, 1, 2);
    diff3 = diff(Fc{3}, 1, 3);
    %Developer Note: Instead of f = diff1 + diff2 + diff3; which calls the 
    % constructor two times, we use the following trick to call it just
    % once. See CHEBFUN3/PLUS for more details:
    vscales = [vscale(diff1) + vscale(diff2), vscale(diff3)];
    
    m = 51; % size of sampling grid
    LVals = sample(diff1, m, m, m) + sample(diff2, m, m, m) + ...
        sample(diff3, m, m, m);
    LVscale = max(abs(LVals(:)));
    kappa = sum(vscales)/LVscale;
    pref = chebfunpref().cheb3Prefs;
    eps = pref.chebfun3eps;
    tol = eps*kappa;
    f = chebfun3(@(x,y,z) feval(diff1, x, y, z) + feval(diff2, x, y, z) + ...
        feval(diff3, x, y, z) , Fc{1}.domain, 'eps', tol);
    
else                      
     error('CHEBFUN:CHEBFUN3V:divergence:notSupported', ...
        'three inputs are needed.')
end

end