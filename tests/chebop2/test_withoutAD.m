function pass = test_withoutAD( pref )
% Try the solvers when the user supplies the AD information and low rank
% representation. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 100*pref.cheb2Prefs.chebfun2eps;


state = warning;
warning('off','all')

% Try Helmoltz:
N = chebop2;
N.domain = [-1 1 -1 1];
N.U = [500 1;0 0;1 0];
N.S = [1 0; 0 1]; 
N.V = [1 500;0 0;0 1];
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
N.xorder = 2; 
N.yorder = 2; 
u = N \ 0; 

N = chebop2(@(u) diffx(u,2) + diffy(u,2) + 1000*u);
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1;  
exact = N \ 0; 

x = linspace(-1,1,1000); 
[xx,yy] = meshgrid(x); 
A = abs(u(xx,yy) - exact(xx,yy));
pass(1) = ( max( A(:) ) < 20000*tol );

% Try variable coefficients:
x = chebfun(@(x) x);
exact = chebfun2(@(x,y) 3*y.^2+x.^3); 
N = chebop2;
N.domain = [-1 1 -1 1];
N.U = {0 1;0 0;1 0};
N.S = [1 0; 0 1]; 
N.V = {-x 0; 0 0; 0 01};
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.ubc = exact(:,1); N.dbc = exact(:,-1); 
N.xorder = 2; 
N.yorder = 2; 
u = N \ 0; 

x = linspace(-1,1,1000); 
[xx,yy] = meshgrid(x); 
A = abs(u(xx,yy) - exact(xx,yy));
pass(2) = ( max( A(:) ) < tol ); 

warning(state);
end