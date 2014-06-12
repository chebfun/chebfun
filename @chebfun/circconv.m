function h = circconv(f,g)
%CIRCCONV   Circular convolution of a Fourier-based chebfun on its interval [a,b].
%   S = CIRCCONV(F,G) is the circular convolution from a to b of F and G.
%   
% See also CONV.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

% No support for quasimatrices:
if ( numColumns(f) > 1 || numColumns(g) > 1 )
    error('CHEBFUN:circconv:quasi', 'No support for array-valued CHEBFUN objects.');
end

% Check transpose state:
if ( xor(f(1).isTransposed, g(1).isTransposed) )
    error('CHEBFUN:circconv:transposed', 'CHEBFUN dimensions do not agree.');
end

if ~isa(f.funs{1}.onefun,'fourtech') && ~isa(g.funs{1}.onefun,'fourtech')
    error('CHEBFUN:circconv:NotAvailable','Circular convolutions only possible for Fourier-based chebfuns.');
end

% Extract the domain:
[a, b] = domain(f);
[c, d] = domain(g);

if a ~=c &&  b ~= d
    error('CHEBFUN:circconv:domain','Domains of f and g must match');
end

% Call BNDFUN circular convolution.
h = chebfun(circconv(f.funs{1}, g.funs{1}),domain(f));

end