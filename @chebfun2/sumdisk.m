function [Integral] = sumdisk(f)
%SUM2DISK   Double integral of a CHEBFUN2 over the unit disk.
%   I = SUM2DISK(F) returns the double integral of a CHEBFUN2 over the unit
%   disk. Chebyshev polynomial approximation formula of maximal degree N is
%   used
%   
% See also INTEGRAL2, INTEGRAL, QUAD2D, SUM2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note for developers: This method supposes that the CHEBFUN2 is based on 
% Chebyshev polynomials.


ChebCoeff_Mat = chebcoeffs2(f);
ChebCoeff_Mat_RowNum = size(ChebCoeff_Mat,1);
ChebCoeff_Mat_ColNum = size(ChebCoeff_Mat,2);

Integral = 0;

for I = 1:ChebCoeff_Mat_RowNum
   for J = 1:ChebCoeff_Mat_ColNum
       Integral = Integral + ChebCoeff_Mat(I,J)*TixTjy_DiscIntegral_Fn(I-1,J-1);
   end
end


end

%%
% This local function computes the integral of T_i(x)T_j(y) over the unit
% disk.

function [ DiscIntegral ] = TixTjy_DiscIntegral_Fn(I,J)
%TixTjy_DiscIntegral_Fn Computes the integral of product of Chebyshev
%polynomials T_i(x)T_j(y) in the unit disc. Slevinsky formula is used
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


