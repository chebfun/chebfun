function [Integral] = sumdisk(f)
%SUM2DISK   Double integral of a CHEBFUN2 over the unit disk.
%   I = SUM2DISK(F) returns the double integral of a CHEBFUN2 over the unit
%   disk. The integral is evaluated using a truncated tensor product 2D Chebyshev
%   polynomial  that resolve the function to
%   machine precision
%
% See also INTEGRAL2, INTEGRAL, QUAD2D, SUM2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note for developers: This method supposes that the CHEBFUN2 is based on
% Chebyshev polynomials.

% Empty check:
if ( isempty( f ) ) 
    Integral = 0;
    return
end

% Check f is a chebfun2 or a trigfun2

colsTechf = get(f.cols.funs{1}, 'tech');

if ( isequal(colsTechf,@trigtech) )
    f = chebfun2(@(x,y) f(x,y));            % slow way to convert a chebfun2 to a trigfun 2
% else ( isequal(colsTechf,@trigtech) )     % No need to check
end

% will add the arbitrary domain later


ChebCoeffMat = chebcoeffs2(f);                  % Matrix of Chebyshev coefficients
ChebCoeffMatRowNum = size(ChebCoeffMat,1);      
ChebCoeffMatColNum = size(ChebCoeffMat,2);

Integral = 0;

% Loop through the entries of the matrix of Chebyshev coefficients
for I = 1:ChebCoeffMatRowNum
    for J = 1:ChebCoeffMatColNum
        Integral = Integral + ChebCoeffMat(I,J)*GetTixTjyDiscIntegral(I-1,J-1);
        % use the local function GETIXTJYDISCINTEGRAL() to compute the
        % integral of T_i(x)T_j(y) over the unit disk.
    end
end


end

%%
% This local function computes the integral of T_i(x)T_j(y) over the unit
% disk.

function [ DiscIntegral ] = GetTixTjyDiscIntegral(I,J)
%TixTjy_DiscIntegral_Fn   Computes the integral of product of Chebyshev
%   polynomials T_i(x)T_j(y) in the unit disc. Slevinsky formula is used
%   here.
%
%   Detailed explanation goes here


if mod(I,2) == 1
    DiscIntegral = 0;
elseif mod(J,2) == 1
    DiscIntegral = 0;
else
    
    if J == 0
        if I == 0
            DiscIntegral = pi;
        elseif I == 2
            DiscIntegral = -pi/2;
        else
            DiscIntegral = 0;
        end
        
    end
    
    
    if J>0
        if I == J-2 && J~=2
            DiscIntegral = (pi/2)*((-1)^(J/2))/(2*(J-1));
        elseif I == J
            DiscIntegral = (pi/2)*(((-1)^(J/2))/(2*(J+1)) - ((-1)^(J/2))/(2*(J-1)));
        elseif I == J+2
            DiscIntegral = -(pi/2)*((-1)^(J/2))/(2*(J+1));
        elseif I == 0 && J == 2
            DiscIntegral = -(pi/2);
        else
            DiscIntegral = 0;
        end
    end
    
    
end

end


