function g = permute(f, dims)
%PERMUTE   for CHEBFUN3 objects. This is a generalization of the transpose 
%   operation for bivariate functions.
%
%   PERMUTE(F, DIMS) permutes the CHEBUFN3 object F in the order specified 
%   by the row vector DIMS. For example, if f is a CHEBFUN3 representation 
%   of f(x, y, z), then
%   permute(f, [1 3 2]) 
%   is a CHEBFUN3 representation of f(x, z, y).

g = f;
dom = f.domain;
if ( dims == [1 2 3] )
    return
elseif ( dims == [1 3 2] )
    temp = g.rows;
    g.rows = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = [dom(1) dom(2) dom(5) dom(6) ...
        dom(3) dom(4)];
    
elseif ( dims == [2 1 3] )
    temp = g.cols;
    g.cols = g.rows;
    g.rows = temp;
    g.core = permute(g.core, dims);
    g.domain = [dom(3) dom(4) dom(1) dom(2) ...
        dom(5) dom(6)];
    
elseif ( dims == [2 3 1] )
    temp = g.cols;
    g.cols = g.rows;
    g.rows = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = [dom(3) dom(4) dom(5) dom(6) ...
        dom(1) dom(2)];
    
elseif ( dims == [3 1 2] )
    temp = g.rows;
    g.rows = g.cols;
    g.cols = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = [dom(5) dom(6) dom(1) dom(2) ...
        dom(3) dom(4)];
    
elseif ( dims == [3 2 1] )
    temp = g.cols;
    g.cols = g.tubes;
    g.tubes = temp;
    g.core = permute(g.core, dims);
    g.domain = [dom(5) dom(6) dom(3) dom(4) ...
        dom(1) dom(2)];
end