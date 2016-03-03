function pass = test_diff( pref )
% Check the diff command in Chebfun2

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6 * pref.cheb2Prefs.chebfun2eps;
j = 1; 

% Battery:  functions, df/dx, and df/dy
f = {@(x,y) x, @(x,y) cos(x).*exp(y), @(x,y) cos(x.*y), @(x,y) x.^2 + x.*y.^2};

fx = {@(x,y) 1+0*x, @(x,y) -sin(x).*exp(y), @(x,y) -y.*sin(x.*y), @(x,y) 2*x+y.^2};

fy = {@(x,y) 0*x, @(x,y) cos(x).*exp(y),@(x,y) -x.*sin(x.*y), @(x,y) 2*x.*y};

% Domains
D = [-1 1 -1 1;... 
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

err = norm( diff(g) - diff(g, 1, 1) );
pass(j) = ( err < tol); j = j + 1;
err = norm( diff(g,1) - diff(g,1,1) );
pass(j) = ( err < tol); j = j + 1;
err = norm( diff(g,2) - diff(g,2,1) );
pass(j) = ( err < 100*tol); j = j + 1;
err = norm( diff(g,2,2) - diff(diff(g,1,2),1,2) );
pass(j) = ( err < 10*tol); j = j + 1;

end