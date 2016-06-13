function I = sumdisk(f)
%SUMDISK   Double integral of a CHEBFUN2 over the unit disk.
%   I = SUMDISK(F) returns the double integral of the CHEBFUN2 F over the
%   unit disk if F is defined on the unit square; otherwise the result is
%   scaled appropriately for a different square or rectangle.  The integral
%   is evaluated by using formulas based on the bivariate Chebyshev 
%   expansion of F.
%
% See also INTEGRAL2, INTEGRAL, QUAD2D, SUM2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) ) 
    I = 0;
    return
end

% If f is a trigfun2, convert it to a chebfun2
Domainf = f.domain;
colsTechf = get(f.cols.funs{1}, 'tech');

if ( isequal(colsTechf,@trigtech) )
    f = chebfun2(@(x,y) f(x,y), Domainf, 'vectorize');         % not as clean as we would like
elseif ( isequal(colsTechf,@chebtech2) ~= 1 )
    error('Error: The input argument is not chebfun2 or trigfun2 object')
end

% Compute the integral

coeffs = chebcoeffs2(f);                 % matrix of Chebyshev coefficients
[nRow, nCol] = size(coeffs);


if ( nCol > 1 && nRow > 1 )
    
    % Extract 0 diagonal of the coeff matrix
    % Note that the command 
    %diag behaves differently when the its input argument is a row/column vector
    CoeffsDiag0 = diag(coeffs,0);               
    
    % Extract the even entries
    CoeffsDiag0 = CoeffsDiag0(1:2:end);
    Diag0Length = length(CoeffsDiag0);
    
    % Compute the integral of T_i(x)T_j(y) over the unit disk for the
    % appropriate i,j
    
    kVec = 2*(0:(Diag0Length-1))';
    Diag0Int = (pi*(-1).^(kVec/2))./(2-2*kVec.^2);
    Diag0Int(1,1) = pi;
    
else
    
    CoeffsDiag0 = coeffs(1,1);
    Diag0Int = pi;
    
end

if ( nCol > 2 )
    
    if ( nRow > 1 )
    CoeffsDiag2 = diag(coeffs,2);
    CoeffsDiag2 = CoeffsDiag2(1:2:end);
    Diag2Length = length(CoeffsDiag2);
    
    kVec = 2*(0:(Diag2Length - 1))';
    Diag2Int = (pi*(-1).^(1+kVec/2))./(4*kVec + 4);
    Diag2Int(1,1) = -pi/2;
    else
        CoeffsDiag2 = coeffs(1,3);
        Diag2Int = -pi/2;
    end
    
else
    CoeffsDiag2 = 0;
    Diag2Int = 0;
end

if ( nRow > 2 )
    
    if ( nCol > 1 )
    CoeffsDiagm2 = diag(coeffs,-2);
    CoeffsDiagm2 = CoeffsDiagm2(1:2:end);
    Diagm2Length = length(CoeffsDiagm2);
    
    kVec = 2*(0:(Diagm2Length - 1))';
    Diagm2Int = (pi*(-1).^(1+kVec/2))./(4*kVec + 4);
    Diagm2Int(1,1) = -pi/2;
    else
        CoeffsDiag2 = coeffs(3,1);
        Diagm2Int = -pi/2;
    end
else
    CoeffsDiagm2 = 0;
    Diagm2Int = 0;
end

% [TODO: what is this?]
% % If nrows and ncols are big enough, you can do the following
% % Of course you'll hto modify this a little in case they're not big enough:
%
% I = pi*C(1,1) - (pi/2)*(C(1,3)+C(3,1);
% for i = 3:2:min(ncols,nrows)
%     I = I + xxx*C(i,i);
%     etc.


I = (Diag0Int')*CoeffsDiag0 + (Diag2Int')*CoeffsDiag2 + (Diagm2Int')*CoeffsDiagm2;


% Rescale the integral for non-default domain


if ( any(Domainf ~= [-1 1 -1 1]) )
    I = I * (Domainf(2) - Domainf(1))*(Domainf(4) - Domainf(3))/4;
end



end

