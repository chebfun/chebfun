function pass = chebop2_eulerTricomi
% Check that Euler--Tricomi equation is working. 
% Alex Townsend, August 2013. 

% First, check that all the particular solutions are being 
% calculated correctly. 

tol = 100*chebfun2pref('eps'); 
j = 1; 

exact = chebfun2(@(x,y) 1+0*x);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1; 

exact = chebfun2(@(x,y) y);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1; 

exact = chebfun2(@(x,y) x);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1; 

exact = chebfun2(@(x,y) x.*y);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1;

exact = chebfun2(@(x,y) 3*y.^2+x.^3);
% exact = chebfun2(@(x,y) y.^2-x.^2);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1;

exact = chebfun2(@(x,y) 3*x.^2+y.^3);
% exact = chebfun2(@(x,y) 3*x.^2+y.^3);
N = chebop2(@(x,y,u) y.*diff(u,2,2) - diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1;

exact = chebfun2(@(x,y) 3*x.^2+y.^3);
% exact = chebfun2(@(x,y) 3*x.^2+y.^3);
N = chebop2(@(x,y,u) -y.*diff(u,2,2) + diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1;

exact = chebfun2(@(x,y) y.^3+x.^3.*y);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1;

exact = chebfun2(@(x,y) 6*x.*y.^2+x.^4);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(j) = ( norm( u - exact ) < tol ); j = j + 1;

end