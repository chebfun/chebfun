function [bcrow, bcvalue] = constructBC(bcArg, bcpos, een, bcn, dom, scl, order)
%CONSTRUCTBC   Discretizes the boundary conditions. 
% 
% INPUTS: 
%   bcArg = linear constraint,
%   bcpos = position of constraint in the other variable,
%   een = discretization size for bcvalue,
%   bcn = discretization size for bcrow,
%   dom = domain of functions that bcArg acts on,
%   scl = length of domain in other variable,
%   order = for periodic constraints only. 
% 
% OUTPUTS: 
%   bcrow = discretized constraint,
%   bcvalue = vector of discretized nonhomogeneous part of the constraint 
%   (satisfying bcrow * X = bcvalue, or X * bcrow' = bcvalue').

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(bcArg, 'chebfun') )
    % Dirichlet conditions are sorted out here: 
    if ( bcpos == -1 )    
        % Dirichlet conditions at x = -1 (or y = -1): 
        bcrow = (-1).^(0:bcn-1).';
    elseif ( bcpos == 1 )
        % Dirichlet conditions at x = 1 (or y = 1): 
        bcrow = ones(bcn, 1);
    else
        % Dirichlet conditions at x = bcpos (or y = bcpos):
        bcrow = cos((0:bcn-1) * acos(bcpos));
    end
    bcvalue = resize(bcArg.coeffs(:), een);
    
elseif ( isa( bcArg, 'function_handle' ) )
    % More general conditions are sorted out here: 
    if ( nargin(bcArg) == 1 )
        % This allows the user to do N.lbc = @(u) u - 1, for example.
        bcArg = @(x,u) u - bcArg(x); 
    end
    if ( nargin(bcArg) == 2 )
        % Assume we have a boundary condition is of the form
        % @(x,u) a*u +b*diff(u) + c*diff(u,2) + ...  + f(x), 
        % where a, b, c, ... are constants.
        
        % Evaluate at a chebfun:
        x = chebfun(@(x) x, dom);
        f = bcArg(x, 0*x);
        
        % Get sizes, depending how the user creates the bc: 
        if ( isa(f, 'chebmatrix') ) 
            nf = size(f, 1); 
        elseif ( isa(f, 'chebfun') ) 
            nf = size(f, 2); 
            f = chebmatrix(mat2cell(f)); 
        end 
        
        bcvalue = zeros(een , nf);
        % Construct f(x):  
        for jj = 1:nf
            g = f{jj};
            % bcvalue = -f as it's going in the RHS: 
            cg = g.coeffs(:);
            bcvalue(:,jj) = -resize(cg, een);
        end
        
        % Now go find the constants in the boundary conditions: 
        L = linearize(chebop(bcArg, dom), [], [], 0, 0);
        p = chebop2.recoverCoeffs(L);
        
        % Set up the boundary rows that will impose the linear constraints: 
        if ( iscell(p) )
            bcrow = zeros(bcn, length(p));
            for jj = 1 : length(p)
                pp = p{jj};
                for kk = 1 : size(pp, 2)
                    c = feval(chebfun(pp(:,kk)), bcpos);
                    if ( abs(c) > 0 )
                        % Each derivative must be scaled by a factor
                        % depends on the other variable's domain. 
                        dx_scaling = abs(2/diff(scl)).^(kk-1);
                        val =  dx_scaling*chebValues(kk-1, bcn, bcpos);
                        bcrow(:, jj) = bcrow(:, jj) + c*val;
                    end
                end
            end
        else
            % If it's simple like Dirichlet, quickly set things up: 
            bcrow = zeros(bcn, 1);
            c = feval(p, bcpos);
            if ( abs(c) > 0 )
                val = chebValues(0, bcn, bcpos);
                bcrow = bcrow + c*val;
            end
        end
    end
    
elseif ( isa(bcArg, 'char') )
    % If the boundary conditions are 'periodic' then try and setup the
    % right bcrows. 
    if ( strcmpi(bcArg, 'periodic') )
        bcrow = zeros(bcn, order);
        bcvalue = zeros(een, order);
        for jj = 1:order 
           bcrow(:,jj) = chebValues(jj-1, bcn, 1) - chebValues(jj-1, bcn, -1);
        end
    else
        error('CHEBFUN:CHEBOP2:constructBC:word', ...
            'Unrecognised boundary condition string.');
    end
else
    error('CHEBFUN:CHEBOP2:constructBC:type', ...
        'Unrecognised boundary condition syntax.');
end

end

function v = resize(v, n)
% Make v of length exact n, by either adding trailing zeros or truncating
% the vector. Return as a column.
if ( length(v) < n )
    v(length(v)+1:n) = 0;
else
    v = v(1:n);
end
v = v(:);
end

function val = chebValues(k, n, x)
%CHEBVALUES   Return the values of Chebyshev {T0^(k)(x),..Tn^(k)(x)}, x being 
% 1 or -1. 
if ( k == 0 )
    val = x.^((0:n-1).');
else
    [ll, kk] = meshgrid((0:n-1), (0:k-1));
    val = (x).^((1:n).').*prod((ll.^2 - kk.^2)./(2*kk+1), 1).';
end

end
