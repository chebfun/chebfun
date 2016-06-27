function f = plus(f, g)
%PLUS   Addition of two CHEBFUN3T obejects.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % CHEBFUN3T + [] = []
    f = [];
    
elseif ( isa(g, 'double') ) % CHEBFUN3T + double
    f.coeffs(1,1,1) = f.coeffs(1,1,1) + g;
    f.vscale = abs(f.vscale + g);
    
elseif ( isa(f, 'double') ) % double + CHEBFUN3T
    % Switch argument order and call PLUS again:
    f = plus(g, f);
    
elseif ( isa(f, 'chebfun3t') && isa(g, 'chebfun3t') )  % CHEBFUN3T + CHEBFUN3T
    if (f.domain == g.domain)
        % Make both CHEBFUN3T objects have the same length:
        fCoeffs = f.coeffs;
        gCoeffs = g.coeffs;
        [nxF, nyF, nzF] = size(fCoeffs);
        [nxG, nyG, nzG] = size(gCoeffs);
        fCoeffsNew = zeros(max(nxF,nxG), max(nyF,nyG), max(nzF,nzG));
        gCoeffsNew = zeros(max(nxF,nxG), max(nyF,nyG), max(nzF,nzG));
        fCoeffsNew(1:nxF, 1:nyF, 1:nzF) = fCoeffs;
        gCoeffsNew(1:nxG, 1:nyG, 1:nzG) = gCoeffs;
        
        % Update coefficients:
        f.coeffs = fCoeffsNew + gCoeffsNew;
        f.vscale = vertscale(f);
    else
        error('CHEBFUN3T:plus: Inputs are not defined on the same domain.');
    end
    
end
end