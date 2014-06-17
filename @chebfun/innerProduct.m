function out = innerProduct(f, g)
%INNERPRODUCT   Compute the inner product of two CHEBFUN objects.
%   INNERPRODUCT(F, G) returns the L2 inner product of the two CHEBFUN objects F
%   and G (conjugate linear in F).
%
%   If F and/or G are array-valued CHEBFUN objects or quasimatrices, then the
%   result is a matrix whose i,j entry is the inner product of the ith column of
%   F with the jth column of G.
%
%   If either F or G is a numeric array, it is cast to a CHEBFUN on the domain
%   of the other argument. The inner product of the resulting CHEBFUN and the
%   other input argument is then computed.
%
% See also NORM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Cast numerical input to a CHEBFUN
if ( isnumeric(f) )
    f = chebfun(f, domain(g));
elseif ( isnumeric(g) )
    g = chebfun(g, domain(f));
end

% If either of the functions is defined on an infinite domain, we need to cast 
% to quasimatrices, since UNBNDFUN INNERPRODUCT does not support array-valued 
% inputs.
if ( any(isinf(domain(f))) || any(isinf(domain(g))) )
    f = cheb2quasi(f);
    g = cheb2quasi(g);
end

numColsF = numColumns(f);
numColsG = numColumns(g);

% Initialise the output:
out = zeros(numColsF, numColsG);

if ( numel(f) == 1 && numel(g) == 1 )
    % Array-valued CHEBFUN case:
    
    % Overlap the CHEBFUN objects:
    [f, g] = overlap(f, g);
    
    % Loop over the FUNs:
    for k = 1:numel(f.funs)
        out = out + innerProduct(f.funs{k}, g.funs{k});
    end 
    
else
    % QUASIMATRIX case:
    
    % Convert to a cell array:
    f = mat2cell(f);
    g = mat2cell(g);
    % Loop over the columns:
    for j = 1:numColsF
        for k = 1:numColsG
            out(j, k) = innerProduct(f{j}, g{k});
        end
    end 

end

end
