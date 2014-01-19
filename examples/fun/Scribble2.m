%% A scribble for Chebfun2 
% Nick Hale and Alex Townsend, 30th August 2013

%%
% (Chebfun2 Example fun/Scribble2.m)
% [Tags: #gift, #SCRIBBLE, #fun, #birthday]

%% Scribble in Chebfun
% The command scribble in Chebfun is one of our favourites. It takes a string
% containing a word or sentence and represents it as a complex-valued piecewise
% smooth chebfun. For instance, it is Nick Trefethen's birthday on Friday so we
% can send him a birthday message:

s = scribble('Happy birthday LNT!');
plot(s,'linewidth',2), axis off, axis equal

%% 
% The scribble command was introduced purely for fun, but over time we found it 
% an entertaining way to illustrate complex variables and to show the ease of 
% piecewise approximation in Chebfun, see [2,3].

%% Scribble in Chebfun2 
% In this Example for Nick's birthday we have developed a scribble for Chebfun2 
% that can be used to write a word or sentence. This time we represent the string 
% by a low-rank function using Chebfun2. For instance, 

s = scribble2('Happy birthday LNT!'); 
contour(s,.3:.05:.85), axis off

%% 
% Thus we find out that the string 'Happy birthday LNT!' can be approximated by 
% a rank 27 function:

rank(s)

%%
% If scribble2 is given no output variable then it plays a movie. The movie
% shows the approximation being constructed by approximations of rank 1, 2, 3,
% and so on. The movie stops when the message has been approximated to machine 
% precision. 

scribble2('Happy birthday LNT!')

%% 
% Historically, in October 2011, this type of movie was one of the first 
% experiments on the Chebfun2 system [1]. Now, as of August 2013, you can repeat 
% the experiment with one line of MATLAB.  

%% References
%% 
% [1] A. Townsend, Hello World, Chebfun2 Example,
% http://www.chebfun.org/examples/fun/html/HelloWorld.shtml, March 2013.
%%
% [2] L. N. Trefethen, Writing a message in 3D, Chebfun Example,
% http://www.chebfun.org/examples/fun/html/Writing3D.shtml, Nov 2010.
%% 
% [3] L. N. Trefethen, Birthday cards and analytic functions, Chebfun Example, 
% http://www.chebfun.org/examples/fun/html/Birthday.shtml, Sept 2010.
