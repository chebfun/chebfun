close all, clear all, clc

d = [-1 0 1];                     % function domain
I = operatorBlock.eye(d);         % identity
D = operatorBlock.diff(d);        % differentiation
x = chebfun(@(x) x, d);           % the variable x on d
E = functionalBlock.eval(d);      % evaluation functional generator
z = functionalBlock.zero(d);

A = [ D^3 I ; 
      D^2 D];
  
A = linop(A);
A = addbc( A, [E(-1), z], 1 );
A = addbc( A, [E(1), z], 0 ); 
A = addbc( A, [z, E(1)], 1 ); 
A = addbc( A, [z, E(-1)], 0 ); 

f = [sin(x) ; cos(x)];
u = A\f;

u = chebfun(u);

figure(1)
plot(u)  

figure(2)
chebpolyplot(u), shg
