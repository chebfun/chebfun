%DISKHARM
%Analogous to sphharm, diskharm(L,m,type)returns is the Lth order, mth eigenfunction (e.g. (J_L(w_m r)*cos(L*theta) ) 
%for laplacian with
%BC specified by "type": 'D' = dirichlet, 'N'=Neumann; default is dirichlet
function Y = diskharm(L,m, type)

if ( nargin < 3 || isempty(type))
    type='D';
end

dom = [-pi pi 0 1];
    Y=diskfun(@(theta, r) mydiskharm(theta, r, L, m, type), 1, dom);
   
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
       end
       
       if type=='N'
           if (L==0 && m==1)
               z = 1+r*0;  %case where eigenvalue is zero; constant mode
           else
            jzero = roots(chebfun(@(x) L*besselj(L,x)-x.^(0^(0^L)).*besselj(L+1, x), [L (m+L/2)*pi])); %the x-term is not present when L=0
            jzero=jzero(m); %when L=0 this is still the correct index since 0 is always the first root and the first counted eigenvalue
            z=besselj(L, r*jzero);
           end
       end
       %choose sin or cos based on L value
       pos = abs(max(0,lsign));  
       Z=((pos)*cos(L*theta) +(1-pos)*sin(L*theta)).*z;
       
       Z=reshape(Z,sz1);
    end  

end