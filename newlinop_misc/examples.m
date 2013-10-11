%%
clc
clear classes

%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros;
x = chebfun('x', dom);
u = chebfun('x.^2', dom);
U = linop.diag(u);   


%% Collocation discretizations
matrix(I,5)
%%
matrix(D,5)
%%
matrix(U,5)
%%
% shorthand notation: U*matrix applies the discretization of U
%v = (1:5)';
%v = [ v flipud(v) ];
%abs( U*v - matrix(U,5)*v)

%% Operator instantiations
blockOp(I)
%%
dop = blockOp(D)
TwoX = dop(x.*x)   % should be the chebfun 2*x
TwoXAgain = D*(x.*x)
ShouldBeZero = norm( 2*x - TwoX )
%%
uop = blockOp(U)
ShouldBeZero = norm( uop(x+1) - x.*x.*(x+1) )   % should be zero
%%
% shorthand notation, U*chebfun
ShouldBeZero = norm( uop(x+1) - U*(x+1) )

%% Coefficient realizations
% Quasimatrix [ 1 ]
Ic = coeff(I)
%%
% Quasimatrix [ 1 0 ]
Dc = coeff(D)
%%
% Quasimatrix [ x.^2 ]
Uc = coeff(U)

%% operator arithmetic
A = D*(U*D) + I
B = D + I;
C = A^3;

%%
% collocation discretization
L = matrix(A,5)
%%
% check it
ShouldBeZero = L - ( diffmat(5)*(diag(x(chebpts(5)).^2)*diffmat(5)) + eye(5) )

%%
% operator
f = op(A)
%%
% check it
ShouldBeZero = norm(f(exp(x)) - (diff(x.*x.*diff(exp(x))) + exp(x)))

%%
% coefficient format
f = coeff(A)
%%
% check it
% ShouldBeZero = norm( f - [x.*x 2*x 1] )

%% chebmatrix
A = chebmatrix( {I,Z;D,U} )

%% 
% realization
matrix(A,4)

%% more complicated chebmatrix
data = { I, x, -I; 
    linop.sum, 5, linop.evalAt(dom,'right');
    D, chebfun(1,dom), U };
A = chebmatrix( data )

%%
% visualization
spy(A)

%%
% realization
matrix(A,4)

%%
% application to an appropriate chebmatrix
vv = { sin(x); pi; cos(x) };
v = chebmatrix( vv );
Av = A*v

%%
% check A*v result
ShouldBeZero = norm( v{1}+pi*x-v{3} - Av{1} )
ShouldBeZero = norm( sum(v{1})+5*v{2}+feval(v{3},dom(2)) - Av{2} )
ShouldBeZero = norm( diff(v{1})+pi+(x.*x).*v{3} - Av{3} )

%% Solve a linear system 
L = chebmatrix( { D, -I; I, D } );
f = chebmatrix( { x; 0*x } );

%%
E = linop.eval(dom);
B1 = [E('left'), -E('right')];
B2 = [linop.sum(dom), E('left')];
L = bc(L,B1,0);
L = bc(L,B2,1);

%%
u = L\f;

%%
plot(u{1},'b'); hold on
plot(u{2},'r'); hold off, shg

%%
% check the ODEs
ShouldBeZero = norm( diff(u{1})-u{2} - f{1} )
ShouldBeZero = norm( u{1} + diff(u{2}) )

%%
% check the BCs
v = u{2};  u = u{1};
ShouldBeZero = abs( u(-2)-v(2) )
ShouldBeZero = abs( sum(u)+v(-2) - 1)

%% Auto-diff demonstration
% Function of two variables u, v. 
u = x.*x;    ujac = chebmatrix({I,Z});
v = cos(x);  vjac = chebmatrix({Z,I});
f1 = u.*u;     J1 = chebmatrix({2*linop.diag(u)}) * ujac;
f2 = f1.*v;    
A1 = chebmatrix({linop.diag(f1)});  A2 = chebmatrix({linop.diag(v)});
J2 = A1*vjac + A2*J1;
f3 = exp(f2);  J3 = chebmatrix({linop.diag(f3)})*J2;
df3du = J3{1};
df3dv = J3{2};

%%
% check via small perturbation
delta = sin(exp(x));  delta = (1e-6/norm(delta))*delta;
difference1 = (exp((u+delta).*(u+delta).*v) - f3);
norm( difference1 - df3du*delta )
difference2 = (exp(u.^2.*(v+delta)) - f3);
norm( difference2 - df3dv*delta )



