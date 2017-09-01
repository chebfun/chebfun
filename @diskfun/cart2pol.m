function f = cart2pol(f, varargin)
%CART2POL    Transforms a DISKFUN F(x,y) in Cartesian coordinates to a 
%   CHEBFUN2 F(th, r) in polar coordinates, with -pi <= th <= pi, and 0 <=
%   r <=1. Working with F in this setting does not ensure that F is smooth
%   at the origin of the disk (i.e. r=0), and can reduce the
%   accuracy of some computations. It may be better to use the 'polar' flag
%   in DISKFUN, or to work directly with the CDR decomposition of F.
% 
%   F = cart2pol(f, 'cdr') returns F(th, r) with -pi <= th <= pi and
%   -1 <= r <= 1. This is equivalent to the CDR decomposition of F. 
% 
% See also CDR.  

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
 
% Use the CDR to build a chebfun2 with trigfuns for rows and chebfuns for
% columns

[c, d, r] = cdr(f);

% Restrict columns down to r in [0 , 1]: 
if nargin < 2
    c = restrict(c, [0 1]); 
end 

% Form a chebfun2
f = c*d*r.';

end