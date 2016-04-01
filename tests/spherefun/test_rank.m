function pass = test_rank( pref )
% Try some pretty functions and ensure k <= min(m,n)

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

j = 1; 

C = 10;
f = @(x,y,z) 1./(1+C*(z.^2 + (x-1).^2 - y.^2).^2);
f = spherefun(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 

C = 100;
f = @(x,y,z) 1./(1+C*((z-1).^2 + x.^2 - y.^2).^2);
f = spherefun(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 

f = @(x,y,z)cos(10*(x.^2+z)).*sin(10*(x+y.^2));
f = spherefun(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 

f = @(x,y,z) tanh(20*x).*tanh(10*y).*cos(50*z.*x.*y+1);
f = spherefun(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end;

end