%% Testing script 

%% Basic constructor 

f = chebfun2(@(x,y) cos(10*x.*y), [-2 2 -2 2]); 
g = chebfun2(@(x,y) cos(10*x.*y), [-2 2 -2 2]); 

%% Basic operations 

h = f + g; 
h = f*g;
h = f.*f
h = -f; 
2*h;

%% Vector-calculus operations 
 
diff(f, 1, 2);
diff(f, 2, 1); 
diff(diff(f, 1, 1), 1, 2) - diff(diff(f, 1, 2),1, 1)


%% Outer-product 

f = chebfun(@(x) sin(x)); 
g = chebfun(@(x) cos(x), [-2 2]); 
h = f * g'; 

%% Vector-calculus
f = chebfun2(@(x,y) cos(10*x.*y), [-2 2 -2 2]); 
F = gradient( f )
f = divergence( F )
f = curl( F ) 
f = cross( F, -F ) 


