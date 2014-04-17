close all, clear all, clc

d = [-1 1];                     % function domain
I = operatorBlock.eye(d);         % identity
D = operatorBlock.diff(d);        % differentiation
x = chebfun(@(x) x, d);           % the variable x on d
E = functionalBlock.eval(d);      % evaluation functional generator
z = functionalBlock.zero(d);
S = functionalBlock.sum(d);

A = [ D^3 I cos(x); 
      D^2 D -x ; 
      S*D E(0) 1];
  
pref = cheboppref;
pref.discretization = @ultraS;
  
A.prefs = pref;
A = linop(A);

A = addbc( A, [E(-1), z,     1], 1 );
A = addbc( A, [E(1),  z,     2], 0 ); 
A = addbc( A, [z,     E(1),  3], 1 ); 
A = addbc( A, [z,     E(-1), 4], 0 ); 

f = [sin(x) ; cos(x) ; pi];
u = A\f;

u = chebfun(u);

figure(1)
plot(u)  

figure(2)
chebpolyplot(u), shg
