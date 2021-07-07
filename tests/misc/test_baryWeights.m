function pass = test_baryWeights(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% NOTE: We only Check the consistency of signs for barycentric weights Greater
% than 0.


% consistency of first Chebyshev points
xpts = 'chebpts'; j = 0; a = 1;

[x, ~, v] = feval(xpts,5,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

% consistency of second Chebyshev points
xpts = 'chebpts';j = 4; a = 2;

[x, ~, v] = feval(xpts,5,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


% consistency of Legendre  points
xpts = 'legpts'; j = 8;

[x, ~, v] = feval(xpts,5);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


[x, ~, v] = feval(xpts,2047);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


% consistency of Legendre-Lobatto points
xpts = 'lobpts'; j = 12;

[x, ~, v] = feval(xpts,5);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


[x, ~, v] = feval(xpts,2047);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


% consistency of Legendre-Radau points
xpts = 'radaupts'; j = 16;

[x, ~, v] = feval(xpts,5);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


% consistency of Hermite points

xpts = 'hermpts'; j = 20;

[x, ~, v] = feval(xpts,5);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


% consistency of Ultraspherical points
xpts = 'ultrapts'; j = 24; a = 1/4;

[x, ~, v] = feval(xpts,5,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 


% consistency of Laguerre points 
xpts = 'lagpts'; j = 28; a = 1/4;

[x, ~, v] = feval(xpts,5,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048,a);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps;


% consistency of Jacobi points 
xpts = 'jacpts'; j =32; a = 1/4; b = 1/2;

[x, ~, v] = feval(xpts,5,a,b);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+1) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,6,a,b);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+2) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2047,a,b);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+3) = norm( sign(v(id))-sign(vb(id)),inf ) < eps; 

[x, ~, v] = feval(xpts,2048,a,b);  
vb = baryWeights(x); id = ( (abs(v) > 0) & (abs(vb) > 0) );
pass(j+4) = norm( sign(v(id))-sign(vb(id)),inf ) < eps;

end
