function pass = test_chebfun2_diff
% Check the diff command in Chebfun2

tol = 1e6 * eps;
j = 1; 

% Battery:  functions, df/dx, and df/dy
f = {@(x,y) x, @(x,y) cos(x).*exp(y), @(x,y) cos(x.*y), @(x,y) x.^2 + x.*y.^2};

fx = {@(x,y) 1+0*x, @(x,y) -sin(x).*exp(y), @(x,y) -y.*sin(x.*y), @(x,y) 2*x+y.^2};

fy = {@(x,y) 0*x, @(x,y) cos(x).*exp(y),@(x,y) -x.*sin(x.*y), @(x,y) 2*x.*y};

% Domains
D = [-1 1 -1 1;... 
     -2 2 -2 2;...
     -3 1 -1 2];
     
% loop over the various options. 
for jj = 1 : length( f )
    for r = 1 : size(D, 1)
        g = chebfun2( f{jj}, D(r,:) ); 
        gx = chebfun2( fx{jj}, D(r,:) );
        gy = chebfun2( fy{jj}, D(r,:) );

        dx = diff(g, 1, 2); 
        dy = diff(g, 1, 1);

        pass(j) = ( norm( dx - gx ) < tol); j=j+1;
        pass(j) = ( norm( dy - gy ) < tol); j=j+1;
        
    end
end

% Check that different syntax for diff. 

jj = 3;  
g = chebfun2( f{jj}, D(r,:) ); 

pass(j) = ( norm( diff(g) - diff(g, 1, 1) ) < tol); j=j+1;
pass(j) = ( norm( diff(g,1) - diff(g,1,1) ) < tol); j=j+1;
pass(j) = ( norm( diff(g,2) - diff(g,2,1) ) < 10*tol); j=j+1;
pass(j) = ( norm( diff(g,2,2) - diff(diff(g,1,2),1,2) ) < 10*tol); j=j+1;

end