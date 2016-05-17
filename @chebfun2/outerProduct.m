function h = outerProduct(f, g)
%OUTERPRODUCT    The outer product of two CHEBFUN objects. 
%
%   H = OUTERPRODUCT(F, G) returns the CHEBFUN2 representing H(x,y) = F(y)G(x),
%   where F and G are two CHEBFUN objects.
%
%   This command is for internal use only. Users are expected to use  f*g' or
%   f*g.'

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This command could use a compression-like algorithm, but instead we will
% just concatenate the rows and columns.

% Empty check: 
if ( isempty( f ) ) 
    h = chebfun2(); 
    return
end

% Empty check: 
if ( isempty( g ) ) 
    h = chebfun2(); 
    return
end

% Check they are both chebfun objects 
if ( ~isa(f, 'chebfun') || ~isa(g, 'chebfun') ) 
   error('CHEBFUN:CHEBFUN2:outerProduct:badInputs', ...
       'Outer product must involve two chebfun objects.');
end

% Extract out domains:
fdom = domain(f);
gdom = domain(g);
dom = [gdom, fdom]; 

h = chebfun2();
% Form outerproduct: 
if ( size(f, 2) == size(g, 1) )
    h.cols = f; 
    h.rows = g.'; 
    h.pivotValues = ones(size(f, 2), 1);
    h.domain = dom; 
else
    error('CHEBFUN:CHEBFUN2:outerProduct:sizes', ...
        'Sizes not consistent for outer product.');
end

end
