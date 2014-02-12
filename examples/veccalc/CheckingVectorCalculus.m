%% Checking vector calculus
% Alex Townsend, 4th March 2013

%% 
% (Chebfun2 example veccalc/CheckingVectorCalculus.m)
% [Tags: #parallelogram law, #gradient theorem, #divgradcurl, #Chebfun2]

%% Introduction
% Chebfun2 is designed to work with vector valued functions defined on
% rectangles, as well as scalar valued ones. Our convention is to use a lower
% case letter for a scalar function, $f$, and an upper case letter for a vector
% function, $F = \left(f_1,f_2\right)^T$. Vector valued functions can be
% approximated by objects of type chebfun2v. There are two standard ways to make
% a chebfun2v

F = chebfun2v(@(x,y) sin(x), @(x,y) sin(y));     % direct construction

f = chebfun2(@(x,y) sin(x)); g = chebfun2(@(x,y) sin(y)); 
G = [f;g];                         % Concatenation of two chebfun2 objects

%% The parallelogram Law 
% Vector addition, denoted $F + G$, yields another chebfun2v and is computed by
% adding the two scalar components together. It satisfies the parallelogram law,
% which can be verified numerically, as in this example:

F = chebfun2v(@(x,y)cos(x.*y),@(x,y)sin(x.*y)); 
G = chebfun2v(@(x,y)x+y,@(x,y)1+x+y);
abs((2*norm(F)^2 + 2*norm(G)^2) - (norm(F+G)^2 + norm(F-G)^2))


%% The gradient theorem
% The gradient of a chebfun2 is represented by a chebfun2v and is a vector that
% points in the direction of steepest ascent of $f$. The gradient theorem says
% that the integral of $grad(f)$ over a curve only depends on the values of $f$
% at the endpoints of that curve. We can check this numerically by using the
% Chebfun2v command integral. This command computes the line integral of a
% vector valued function. Here we check one example to confirm that if the
% gradient theorem holds:

f = chebfun2(@(x,y)sin(2*x)+x.*y.^2);        % chebfun2 object
F = grad(f);                                 % gradient (chebfun2v)
C = chebfun(@(t) t.*exp(100i*t),[0 pi/10]);  % spiral curve
v = integral(F,C); ends = f(pi/10,0)-f(0,0); % line integral
abs(v-ends)                                  % gradient theorem

%% Curl of the gradient 
% If the chebfun2v $F$ describes a vector velocity field of fluid flow, then
% $curl(F)$ is a scalar function equal to twice the angular speed of a particle
% in the flow at each point. A particle moving in a gradient field has zero
% angular speed and hence, $curl(grad(f))=0$, a well known identity that can
% also be checked numerically.  For example,

norm(curl(grad(f)))

%% More information 
% The code found in this Example can also be found in [1] along with additional
% information about vector calculus in Chebfun2.

%% References
% [1] A. Townsend and L. N. Trefethen, An extension of Chebfun to two 
%     dimensions, submitted (2013). 