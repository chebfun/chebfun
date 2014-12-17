About
=====

Chebfun is an open-source software system for numerical computing with
functions. The mathematical basis of Chebfun is piecewise polynomial
interpolation implemented with what we call “Chebyshev technology”. The
foundations are described, with Chebfun examples, in the book _Approximation
Theory and Approximation Practice_. Chebfun has extensive capabilities for
dealing with linear and nonlinear differential and integral operators, and it
also includes continuous analogues of linear algebra notions like QR and
singular value decomposition. The Chebfun2 extension works with functions of
two variables defined on a rectangle in the x-y plane. To get a sense of the
breadth and power of Chebfun, a great place to start is by looking at our
[Examples][1].


Installation and requirements
=============================

Chebfun is compatible with MATLAB 7.8 (R2009a) and later.

To install, you can either clone the directory with Git or download a .zip
file. Note that a call to `clear classes` is required if you had a previous
version of Chebfun installed.

## Option 1: Download .zip file

Download a .zip of Chebfun from

- https://github.com/chebfun/chebfun/archive/master.zip

After unzipping, you will need to add Chebfun to the MATLAB path. You can do
this either (a) by typing
```
addpath(chebfunroot), savepath
```
where `chebfunroot` is the path to the unzipped directory, (b) by selecting the
`chebfun` directory with the `pathtool` command, or (c) though the File > Set
Path... dialog from the MATLAB menubar.

## Option 2: Clone with Git

To clone the Chebfun repository, first navigate in a terminal to where you
want the repository cloned, then type
```
git clone https://github.com/chebfun/chebfun.git
```
To use Chebfun in MATLAB, you will need to add the `chebfun` directory
to the MATLAB path as above.


Getting started
===============

We recommend taking a look at the [Chebfun Guide][2] and the [Examples
collection][1]. The Guide is an in-depth tour of Chebfun's mathematical
capabilities. The Examples, which number well over one hundred, illustrate
everything from rootfinding to optimization to nonlinear differential
equations and vector calculus. Many users use the Examples as templates for
their own problems.

To get a taste of what computing with Chebfun is like, type
```
x = chebfun('x');
```
and start playing. The variable `x` is a chebfun and can be manipulated in a
way that feels symbolic, although everything Chebfun does is numeric. So try,
for instance:
```
f = sin(12*x).*exp(-x);         % A function on [-1, 1]
g = max(f, 1./(x+2));           % The max of f and 1./(x+2)
plot(g)                         % A function with discontinuous derivative
sum(g)                          % The integral of g
plot(diff(g))                   % The derivative of g
h = g + x - .8;                 % A function with several roots in [-1, 1]
rr = roots(h);                  % Compute the roots of h
plot(h, 'k', rr, h(rr), 'ro')   % Plot h and its roots
```


License
=======

See `LICENSE.txt` for Chebfun's licensing information.



[1]: http://www.chebfun.org/examples/
[2]: http://www.chebfun.org/docs/guide/
