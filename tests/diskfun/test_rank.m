function pass = test_rank()
% Try some pretty functions and ensure k <= min(m,n)

C = 10;
f = @(x,y) 1./(1+C*((x-1).^2 - y.^2).^2);
f = diskfun(f); 
k = length(f); 
[m,n] = length(f); 
if ( k <= min(m,n) )
    pass(1) = 1; 
else
    pass(1) = 0;
end

C = 100;
f = @(x,y) 1./(1+C*( x.^2 - y.^2).^2);
f = diskfun(f); 
k = length(f); 
[m,n] = length(f); 
if ( k <= min(m,n) )
    pass(2) = 1; 
else
    pass(2) = 0;
end

f = @(x,y) cos(10*(x.^2)).*sin(10*(x+y.^2));
f = diskfun(f); 
k = length(f); 
[m,n] = length(f); 
if ( k <= min(m,n) )
    pass(3) = 1; 
else
    pass(3) = 0;
end

f = @(x,y) tanh(20*x).*tanh(10*y).*cos(50*x.*y+1);
f = diskfun(f); 
k = length(f); 
[m,n] = length(f); 
if ( k <= min(m,n) )
    pass(4) = 1; 
else
    pass(4) = 0;
end

end