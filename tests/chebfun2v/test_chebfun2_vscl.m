function pass = chebfun2_vscl
% Check correct vertical scaling. 
% Alex Townsend, March 2013. 

% Laplace operator.
x = chebfun2(@(x,y) x);
try
    f = diff(x,2,1) + diff(x,2,2);
    pass(1)=1;
catch
   pass(1)=0; 
end

% check mdivide
f = x/2; g=x./2;
pass(2) = (abs(f(1,1) - .5)<1e-14);
pass(3) = (abs(g(1,1) - .5)<1e-14);
pass(4) = (abs(f.scl - .5)<1e-14);  % scl getting changed correctly?
pass(5) = (abs(g.scl - .5)<1e-14);

f = 2*x; g = 2.*x;
pass(6) = (abs(f(1,1) - 2)<1e-14);
pass(7) = (abs(g(1,1) - 2)<1e-14);
pass(8) = (abs(f.scl - 2)<1e-14);
pass(9) = (abs(g.scl - 2)<1e-14);

end