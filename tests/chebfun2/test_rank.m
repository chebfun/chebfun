function pass = test_rank( pref )
% Try some pretty functions and ensure k <= min(m,n)

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = pref.eps; 
j = 1; 

C = 100;
f = @(x,y) 1./(1+C*(x.^2 - y.^2).^2);
f = chebfun2(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 

C = 100;
f = @(x,y) 1./(1+C*(.5 - x.^2 - y.^2).^2);
f = chebfun2(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 


C = 1000;
f = @(x,y) 1./(1+C*((x-.5).^2.*(y+.5).^2.*(x+.5).^2.*(y-.5).^2));
f = chebfun2(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 


f = @(x,y)cos(10*(x.^2+y)).*sin(10*(x+y.^2));
f = chebfun2(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 


f = @(x,y) real(airy(5*(x + y.^2))).*real(airy(-5*(x.^2+y.^2)));
f = chebfun2(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end; j = j + 1; 


f = @(x,y) tanh(10*x).*tanh(10*y)./tanh(10).^2 + cos(5*x);
f = chebfun2(f); k = length(f); [m,n]=length(f); 
if k <= min(m,n), pass(j) = 1 ; else pass(j)=0;end;

end