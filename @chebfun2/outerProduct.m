function h = outerProduct(f, g)
%OUTERPRODUCT  The outer product of two chebfun2 objects. 
% 
% This command is for internal use only. Users are expected to use  f*g' or 
% f*g.'
% 
% H = OUTERPRODUCT(F,G) returns the chebfun2 representing H(x,y) = F(y)G(x),
% where F and G are two chebfuns. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


% This command could use a compression-like algorithm, but instead we will just
% concatenate the rows and columns. 

% Empty check: 
if ( isempty( f ) ) 
    h = chebfun2; 
    return
end

% Empty check: 
if ( isempty( g ) ) 
    h = chebfun2; 
    return
end

% Check they are both chebfun objects 
if ( ~isa(f, 'chebfun') || ~isa(g, 'chebfun') ) 
   error('CHEBFUN2:OUTERPRODUCT:INPUTS',...
                            'Outer product must involve two chebfun objects.');
end

% Extract out domains:
fdom = f.domain; 
gdom = g.domain; 
dom = [gdom fdom]; 

h = chebfun2(0, dom); 
% Form outerproduct: 
if ( size(f, 2) == size(g, 1) )
    h.cols = f; 
    h.rows = g.'; 
    h.pivotValues = ones(size(f, 2), 1);
    h.domain = dom; 
else
    error('CHEBFUN2:OUTERPRODUCT:SIZES', 'Sizes not consistent for outerproduct');
end

end