function f = cart2pol(f)
%CART2POL    Transforms a diskfun F(x,y) in Cartesian coordinates to a 
% Chebfun2 F(th, r) in polar coordinates, with -pi <= th <= pi, 
% and 0 <= r <=1. Working with F in this setting does not ensure that the 
% structure associated with the disk is preserved, and can reduce the
% accuracy of some computations. It may be better to use the 'polar' flag in
% Diskfun, or to work directly with the CDR decomposition of F. 
%
% See also CDR.  

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
 
%use CDR to build a chebfun2 with trigfuns for rows and chebfuns for
%columns


[c, d, r] = cdr(f);

%restrict columns down to r in [0 , 1]: 

c = restrict(c, [0 1]); 

%form a chebfun2

f = c*d*r.';

end