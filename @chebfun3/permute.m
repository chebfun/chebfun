function g = permute(f, dims)
%PERMUTE   Permutation of CHEBFUN3 objects. 
%   PERMUTE(F, DIMS) rearranges the dimensions of the CHEBUFN3 object F
%   in the order specified by the row vector DIMS.
%   This is a generalization of the transpose operation for bivariate 
%   functions. For example, if F is a CHEBFUN3 representation of f(x, y, z)
%   then permute(f, [1 3 2]) is a CHEBFUN3 representation of f(x, z, y).
%
% See also PERMUTE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty CHEBFUN3.
if ( isempty(f) )
    g = chebfun3(); 
    return
end

g = f;
dom = f.domain;
if ( dims == [1 2 3] )
    return
elseif ( dims == [1 3 2] )
    temp = g.rows;
    g.rows = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = dom([1 2 5 6 3 4]);
    
elseif ( dims == [2 1 3] )
    temp = g.cols;
    g.cols = g.rows;
    g.rows = temp;
    g.core = permute(g.core, dims);
    g.domain = dom([3 4 1 2 5 6]);
    
elseif ( dims == [2 3 1] )
    temp = g.cols;
    g.cols = g.rows;
    g.rows = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = dom([3 4 5 6 1 2]);
    
elseif ( dims == [3 1 2] )
    temp = g.rows;
    g.rows = g.cols;
    g.cols = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = dom([5 6 1 2 3 4]);
    
elseif ( dims == [3 2 1] )
    temp = g.cols;
    g.cols = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = dom([5 6 3 4 1 2]);
end
