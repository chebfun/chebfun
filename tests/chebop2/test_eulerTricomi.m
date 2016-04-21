function pass = test_eulerTricomi( prefs )
% Check that Euler--Tricomi equation is working. 
% Alex Townsend, August 2013. 

% First, check that all the particular solutions are being 
% calculated correctly. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.techPrefs.chebfuneps; 

exact = chebfun2(@(x,y) 1+0*x);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(1) = ( norm( u - exact ) < tol ); 

exact = chebfun2(@(x,y) y);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(2) = ( norm( u - exact ) < tol );

exact = chebfun2(@(x,y) x);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(3) = ( norm( u - exact ) < tol ); 

exact = chebfun2(@(x,y) x.*y);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(4) = ( norm( u - exact ) < tol );

exact = chebfun2(@(x,y) 3*y.^2+x.^3);
%exact = chebfun2(@(x,y) y.^2-x.^2);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0;
pass(5) = ( norm( u - exact ) < tol ); 

exact = chebfun2(@(x,y) 3*x.^2+y.^3);
% exact = chebfun2(@(x,y) 3*x.^2+y.^3);
N = chebop2(@(x,y,u) y.*diff(u,2,2) - diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(6) = ( norm( u - exact ) < tol );

exact = chebfun2(@(x,y) 3*x.^2+y.^3);
% exact = chebfun2(@(x,y) 3*x.^2+y.^3);
N = chebop2(@(x,y,u) -y.*diff(u,2,2) + diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(7) = ( norm( u - exact ) < tol );

exact = chebfun2(@(x,y) y.^3+x.^3.*y);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(8) = ( norm( u - exact ) < tol );

exact = chebfun2(@(x,y) 6*x.*y.^2+x.^4);
N = chebop2(@(x,y,u) diff(u,2,2) - x.*diff(u,2,1)); 
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
u = N \ 0 ;
pass(9) = ( norm( u - exact ) < tol ); 

end