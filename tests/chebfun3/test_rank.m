function pass = test_rank()
% Try some functions and ensure rank is as expected

C = 1;
ff = @(x,y,z) 1./(1+C*(x.^2 - y.^2 + z.^2).^2);
f = chebfun3(ff); 
r = rank(f); 
m =length(f);
if r <= 30
    pass(1) = 1 ; 
else
    pass(1)=0;
end; 

if m <= 100
    pass(2) = 1 ; 
else
    pass(2)=0;
end;

C = 10;
ff = @(x,y,z) 1./(1+C*(x.^2 - y.^2 + z.^2).^2);
f = chebfun3(ff);
r = rank(f); 
m =length(f);
if r <= 60
    pass(3) = 1 ; 
else
    pass(3)=0;
end; 

if m <= 200
    pass(4) = 1 ; 
else
    pass(4)=0;
end; 

ff = @(x,y,z) real(airy(5*(x + y.^2 + z.^2))) .* ...
    real(airy(-5*(x.^2+y.^2 + z.^2)));
f = chebfun3(ff); 
r = rank(f); 
m =length(f);
if r <= 60
    pass(5) = 1 ; 
else
    pass(5)=0;
end; 

if m <= 200
    pass(6) = 1 ; 
else
    pass(6)=0;
end

end