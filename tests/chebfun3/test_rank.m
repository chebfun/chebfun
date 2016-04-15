function pass = test_rank()
% Try some functions and ensure rank is as expected

j = 1;
C = 1;
ff = @(x,y,z) 1./(1+C*(x.^2 - y.^2 + z.^2).^2);
f = chebfun3(ff); 
r = rank(f); 
m =length(f);
if r <= 30
    pass(j) = 1 ; 
else
    pass(j)=0;
end; 
j = j + 1;
if m <= 100
    pass(j) = 1 ; 
else
    pass(j)=0;
end; 
j = j + 1;

C = 10;
ff = @(x,y,z) 1./(1+C*(x.^2 - y.^2 + z.^2).^2);
f = chebfun3(ff);
r = rank(f); 
m =length(f);
if r <= 60
    pass(j) = 1 ; 
else
    pass(j)=0;
end; 
j = j + 1;
if m <= 200
    pass(j) = 1 ; 
else
    pass(j)=0;
end; 
j = j + 1;

ff = @(x,y,z) real(airy(5*(x + y.^2 + z.^2))) .* ...
    real(airy(-5*(x.^2+y.^2 + z.^2)));
f = chebfun3(ff); 
r = rank(f); 
m =length(f);
if r <= 60
    pass(j) = 1 ; 
else
    pass(j)=0;
end; 
j = j + 1;
if m <= 200
    pass(j) = 1 ; 
else
    pass(j)=0;
end;

end