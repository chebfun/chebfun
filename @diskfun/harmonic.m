function Y = harmonic(L,m, type)
% HARMONIC  Normalized, real-valued, fixed-height cylindrical harmonic.
%
% The cylindrical harmonics are an orthogonal basis of functions in cylindrical 
% coordinates. If the height variable is set as a fixed value, then cylindrical 
% harmonics are the eigenfunctions of Laplace's equation on the disk. 
% For Dirichlet boundary conditions, they are of the form 
% V_L^M(t, r) = A*exp(iLt).*J_L(a_Mr), 
% where J_L is the Lth j-Bessel function, a_M is the Mth positive zero 
% of the function, and A is a normalization factor.
%
% Y = harmonic(L, M) returns the cylindrical harmonic V_L^M(t,r) with Dirichlet
% boundary conditions. Here,
%        -pi <= t <= pi   is the angular coordinate, and
%          0 <= r  <= 1   is the radial coordinate.
%
% Y = harmonic(L, M,'N') returns the cylindrical harmonic V_L^M(t,r) with Neumann
% boundary conditions. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TO DO: This code can be made much faster.


if ( nargin < 3 || isempty(type))
    type='D';
end

dom = [-pi pi 0 1];
    Y=diskfun(@(theta, r) mydiskharm(theta, r, L, m, type), 'polar');
   
%create disk harmonic function 

    function   Z=mydiskharm(theta, r, L, m, type)
        
       sz1=size(theta);
       theta=theta(:);
       r=r(:);
       
       %calculate the mth positive zero of the Lth Bessel function
       
       lsign=sign(L+1); 
       L=abs(L);
       
       %searches for mth zero using lower bound on m=1 and McMahon's asymptotic upper bound
       % see R. C. McCann, "Lower bounds for the zeros of Bessel functions"
       if type=='D'
       jzero = roots(chebfun(@(x) besselj(L,x), [sqrt((3/4)^2*pi^2+L^2) (m+L/2)*pi]));
       jzero=jzero(m);
       z = besselj(L,r*jzero);
       %choose sin or cos based on L value
       pos = abs(max(0,lsign));  
       Z=((pos)*cos(L*theta) +(1-pos)*sin(L*theta)).*z;
       %reshape and normalize
       Z=reshape(Z,sz1);
       k = double(L==0);
       Z = sqrt(2)/ (sqrt((1+k)*pi)*abs(besselj(L+1, jzero)))*Z; 
       end
       
       if type=='N'
           if (L==0 && m==1)
               z = 1+r*0;  %case where eigenvalue is zero; constant mode
              Z = reshape(z, sz1);
              return ; 
           else
            jzero = roots(chebfun(@(x) L*besselj(L,x)-x.^(min(1,L))...            
                .*besselj(L+1, x), [L (m+L/2)*pi])); %the x-term is not present when L=0
            jzero=jzero(m); %when L=0 this is still the correct index 
            z=besselj(L, r*jzero);
           end
           %choose sin or cos based on L value
       pos = abs(max(0,lsign));  
       Z=((pos)*cos(L*theta) +(1-pos)*sin(L*theta)).*z;
       %reshape and normalize
       Z=reshape(Z,sz1);
       k = double(L==0);
       Z = sqrt(2)/ (sqrt( (1-L^2/jzero^2)*(1+k)*pi)*abs(besselj(L, jzero)))*Z; 
       end
        
    end  
    
  
end