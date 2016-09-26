function k = outerProd(f, g, h)
%OUTERPRODUCT    The outer product of three CHEBFUN objects. 
%
% This should be moved to chebfun/outerProduct

%   K = OUTERPRODUCT(F, G, H) returns the CHEBFUN3 representing 
%   H(x,y) = F(x)G(y)H(z),
%   where F, G and H are three CHEBFUN objects.

% TODO: This command could use a compression-like algorithm, but instead 
%   here we are just concatenating cols, rows and tubes.

% Empty check: 
if ( isempty(f) || isempty(g) || isempty(h) ) 
    k = chebfun3(); 
    return
end

% Check they are all chebfun objects 
if ( ~isa(f, 'chebfun') || ~isa(g, 'chebfun') || ~isa(h, 'chebfun')) 
   error('CHEBFUN:CHEBFUN3:outerProd:badInputs', ...
       'Outer product must involve three chebfun objects.');
end

% Extract out information:
fdom = domain(f);
gdom = domain(g);
hdom = domain(h);
dom = [fdom, gdom, hdom]; 

k = chebfun3();
% Form outerproduct: 
if ( size(f, 2) == size(g, 2) && size(f, 2) == size(h, 2) )
    k.cols = f; 
    k.rows = g; 
    k.tubes = h;
    k.core = ones(size(f, 2), size(g, 2), size(h, 2));
    k.domain = dom; 
else
    error('CHEBFUN:CHEBFUN3:outerProd:sizes', ...
        'Sizes not consistent for outer product.');
end

end
