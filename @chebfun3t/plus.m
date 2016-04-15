function f = plus(f,g)
%PLUS    addition of two chebfun3t obejects

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % CHEBFUN3T + [] = []
    
    f = [];
    
elseif ( isa(g, 'double') ) % CHEBFUN3T + double
    f.coeffs(1,1,1) = f.coeffs(1,1,1) + g;
    f.vscale = f.vscale + g;
    %f.epslevel = f.epslevel;
    
elseif ( isa(f, 'double') ) % double + CHEBFUN3T
    % Switch argument order and call PLUS again:
    f = plus(g, f);
    
elseif ( isa(f, 'chebfun3t') && isa(g, 'chebfun3t') )  % CHEBFUN3T + CHEBFUN3T
    if (f.domain == g.domain)
        
        % Make both CHEBFUN3T objects have the same length:
        fcoeffs = f.coeffs;
        gcoeffs = g.coeffs;
        [nxF, nyF, nzF] = size(fcoeffs);
        [nxG, nyG, nzG] = size(gcoeffs);
        fcoeffsNew = zeros(max(nxF,nxG), max(nyF,nyG), max(nzF,nzG));
        gcoeffsNew = zeros(max(nxF,nxG), max(nyF,nyG), max(nzF,nzG));
        fcoeffsNew(1:nxF, 1:nyF, 1:nzF) = fcoeffs;
        gcoeffsNew(1:nxG, 1:nyG, 1:nzG) = gcoeffs;
        
        % Update coefficients:
        f.coeffs = fcoeffsNew + gcoeffsNew;
        f.vscale = f.vscale + g.vscale;
        %f.epslevel = max(f.epslevel , g.epslevel);
else
    error('CHEBFUN3T:plus: Inputs are not defined on the same domain.\n');
end
    
%     f.coeffs = f.coeffs + g.coeffs;
%     % Look for a zero output:
%     tol = max(f.epslevel.*f.vscale, g.epslevel.*g.vscale);
%     absCoeffs = abs(f.coeffs);
%     isz = bsxfun(@lt, absCoeffs, .2*tol); % Are coeffs below .2*el*vs?
%     
%     if ( all(isz(:)) )
%         % Create a zero CHEBTECH:
%         epslevel = max(f.epslevel, g.epslevel);
%         ishappy = f.ishappy && g.ishappy;
%         z = zeros(1, size(f.coeffs, 2));
% 
%         data.vscale = z;
%         data.hscale = f.hscale;
%         f = f.make(z, data);
%         f.epslevel = epslevel;
%         f.ishappy = ishappy;
%     else
%         % Update vscale, epslevel, and ishappy:
%         vscaleNew = getvscl(f); 
%         % See CHEBTECH CLASSDEF file for documentation on this:
%         tmpVscaleNew = vscaleNew;
%         tmpVscaleNew(tmpVscaleNew == 0) = 1;  % Avoid NaNs.
%         epslevelBound = (f.epslevel.*f.vscale + g.epslevel.*g.vscale)./tmpVscaleNew;
%         f.epslevel = updateEpslevel(f, epslevelBound);
%         f.vscale = vscaleNew;

%        f.ishappy = f.ishappy && g.ishappy;
end
end