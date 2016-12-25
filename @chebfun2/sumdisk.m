function I = sumdisk(f)
%SUMDISK   Double integral of a CHEBFUN2 over the unit disk.
%   I = SUMDISK(F) returns the double integral of the CHEBFUN2 F over the
%   unit disk if F is defined on the unit square; otherwise the result is
%   scaled appropriately for a different square or rectangle.  The integral
%   is evaluated by using formulas based on the bivariate Chebyshev 
%   expansion of F.
%
% See also INTEGRAL2, INTEGRAL, QUAD2D, SUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) ) 
    I = 0;
    return
end

% If f is a trigfun2, convert it to a chebfun2
Domainf = f.domain;
colsTechf = get(f.cols.funs{1}, 'tech');

if ( isequal(colsTechf,@trigtech) )
    % not as clean as we would like
    f = chebfun2(@(x,y) feval(f, x, y), Domainf, 'vectorize');
end

% Extract even-index Chebyshev coefficients
coeffs = chebcoeffs2(f);              
coeffs = coeffs(1:2:end, 1:2:end); 
[nRow, nCol] = size(coeffs);

% Entries on the main diagonal
if ( nRow == 1 )
    CoeffsDiag0 = coeffs(1, 1);
    Diag0Int = pi;
    
elseif ( nCol == 1 ) 
    CoeffsDiag0 = coeffs(1,1);
    Diag0Int = pi;
    
else    % nRow > 1 and nCol > 1
    % Extract 0 diagonal of the coeff matrix
    CoeffsDiag0 = diag(coeffs, 0);               
    Diag0Length = length(CoeffsDiag0);
    
    % Compute the integral of T_i(x)T_j(y) over the unit disk for the
    % appropriate i,j
    kVec = 2*(0:(Diag0Length-1))';
    Diag0Int = (pi*(-1).^(kVec/2))./(2-2*kVec.^2);
    Diag0Int(1, 1) = pi;
end

% Entries on the +2 diagonal

if ( nCol == 1 )
    CoeffsDiag2 = 0;
    Diag2Int = 0;
    
else   % nCol > 1
    
    if ( nRow == 1 )
        CoeffsDiag2 = coeffs(1, 2);
        Diag2Int = -pi/2;
    else        % nRow > 1
        CoeffsDiag2 = diag(coeffs, 1);
        Diag2Length = length(CoeffsDiag2);
        kVec = 2*(0:(Diag2Length - 1))';
        Diag2Int = (pi*(-1).^(1+kVec/2))./(4*kVec + 4);
        Diag2Int(1,1) = -pi/2;
    end
    
end

% Entries on the -2 diagonal
if ( nRow == 1 )
    CoeffsDiagm2 = 0;
    Diagm2Int = 0;
    
else   % nRow > 1
    if ( nCol == 1 )
        CoeffsDiag2 = coeffs(3,1);
        Diagm2Int = -pi/2;
    else        % nRow > 1
        CoeffsDiagm2 = diag(coeffs,-1);
        Diagm2Length = length(CoeffsDiagm2);
        kVec = 2*(0:(Diagm2Length - 1))';
        Diagm2Int = (pi*(-1).^(1+kVec/2))./(4*kVec + 4);
        Diagm2Int(1,1) = -pi/2;
    end
end

% Sum up all the terms
I = (Diag0Int')*CoeffsDiag0 + (Diag2Int')*CoeffsDiag2 + ...
    (Diagm2Int')*CoeffsDiagm2;

% Rescale the integral for non-default domain
if ( any(Domainf ~= [-1 1 -1 1]) )
    I = I * (Domainf(2) - Domainf(1))*(Domainf(4) - Domainf(3))/4;
end

end