%% Hello World
% Alex Townsend, 4th March 2013

%% 
% (Chebfun2 Example: Fun/HelloWorld.m)
% [Tags: #Chebfun2, #Rank]

%%
% In any programming language printing "Hello World" to the command line
% is always a first example.  Here we display "Hello" using
% Chebfun2. Can you adapt it to display "Hello World" instead? 

%% 
% Here is a matrix that encodes the word "Hello", from Exercise 9.3 of [1]. 
A=zeros(15,40); 
A(2:9,2:3)=1; A(5:6,4:5)=1;A(2:9,6:7)=1; A(3:10,10:11)=1;
A(3:4,10:15)=1; A(6:7,10:15)=1; A(9:10,10:15)=1; A(4:11,18:19)=1;
A(10:11,18:24)=1; A(5:12,26:27)=1; A(11:12,26:31)=1;
A(6:13,34:35)=1; A(6:13,38:39)=1; A(6:7,36:37)=1; A(12:13,36:37)=1;
spy(A)  % spy plot

%%
% The matrix is size 15 by 40 and hence, of rank at most 15.  In this case
% it is of rank 10 because there are five zero rows:

rank(A)

%% Constructing a chebfun2 from discrete data
% Usually Chebfun2 is passed a function of two variables, but it can also
% deal with discrete data such as a matrix, with syntax such as 
% chebfun2(A). The matrix A, of size m by n, is assumed 
% to contain data values of a function sampled at a m by n Chebyshev 
% tensor grid, and the resulting chebfun2 interpolates A, for example:

f = chebfun2(A);           % chebfun2
X = chebpolyval2(f);       % evaluate on a grid
norm(A - X)                % interpolation error

%% Saying Hello
% We can also pass the Chebfun2 constructor an integer k so that the 
% resulting chebfun2 is of rank exactly k. Here is one way to 
% say "Hello": 
    
m = 200;
x = linspace(-1,1,m);
[xx yy]=meshgrid(x);
[ss tt]=chebpts2(m);

B = flipud(A);             % flip because of matrix indexing
for k = [1 3 5 7 10]
    f = chebfun2(B,k);
    X = f(ss,tt);
    figure
    contour(xx,yy,X,.1:.1:.99), axis off
    title(sprintf('Rank %u',k),'fontsize',16)
    pause(.1)
end

%% 
% References:
%% 
% [1] L. N. Trefethen and D. Bau III, Numerical Linear Algebra, SIAM, 1997.