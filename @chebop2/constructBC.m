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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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
        % Do we have bcs with derivatives in them.
        try
            % This works if the BC operators are vectorized. 
            % Evaluate at xy for AD information: 
            fcell = bcArg(adchebfun2(chebfun2(0, [dom, dom])),...
                           abchebfun2(chebfun2(@(x,y) x.*y, [dom, dom])));
        catch
            % Cannot do it all at once, so do it term by term: 
            gg = adchebfun2(chebfun2(@(x,y) cos(x.*y), [dom, dom]));
            fcell = cell(1, nf); 
            for jj = 1:nf
                fcell{jj} = diff(gg, jj-1);
            end
        end
        
        % Find constants: 
        cc = ones(1,nf);
        for jj = 1:nf
            if ( isa(fcell, 'cell') )
                v = fcell{jj};
            else
                v = fcell;
            end
            cc(jj) = abs(diff(scl)/2).^(length(v.jacobian) - 1);
        end
        
        % Find f(x):  
        for jj = 1:nf
            g = f{jj};
            % bcvalue = -f as it's going in the RHS: 
            bcvalue(:,jj) = -resize(cc(jj)*g.coeffs(:), een);
        end
        
        % Now go find the constants in the boundary conditions: 
        L = linearize(chebop(bcArg, dom));
        p = recoverCoeffs(L);
        
        % Set up the boundary rows that will impose the linear constraints: 
        if ( iscell(p) )
            bcrow = zeros(bcn, length(p));
            for jj = 1 : length(p)
                pp = p{jj};
                for kk = 1 : size(pp, 2)
                    c = feval(chebfun(pp(:,kk)), bcpos);
                    if ( abs(c) > 0 )
                        val = chebValues(kk-1, bcn, bcpos);
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

function p = recoverCoeffs(L)
%RECOVERCOEFFS   Recover coefficient functions of a linear operator.
%   P = RECOVERCOEFFS(L) returns, for a linear operator L, a CHEBFUN
%   quasimatrix P such that
%         Lu = P(:,1)*u + P(:,2)*u' + P(:,3)*u" + ... P(:,M+1)*u^(M),
%   where M is the difforder of the operator. If L is not linear, an error is
%   thrown.
%
%    For a block operator L, i.e., one defining a system of equations
%         Lu = [L_{1,1} L_{1,2} ... L_{1,S}] [ u_1 ]
%              [L_{2,1} L_{1,2} ... L_{1,S}] [ u_2 ]
%              [  ...     ...   ...   ...  ] [ ... ]
%              [L_{R,1} L_{R,2} ... L_{R,S}] [ u_S ],
%    P will be the RxS cell array such that P{J,K} = RECOVERCOEFFS(L_{J,K}).
%
%    [P L] = RECOVERCOEFFS(L) returns also the linop L, which can be useful if
%    the input was a linear chebop.
%
% Example 1:
%  [L x] = chebop(@(x,u) 0.5*diff(u,2) - sin(x).*diff(u) + x.*u);
%  p = recoverCoeffs(L)
%  norm(p - [x -sin(x) 0.5])
%
% Example 2:
%  [L x] = chebop(@(x,u) diff(sin(x).*(diff(cos(x).*u))),[-pi pi]);
%  p = recoverCoeffs(L)
%  norm(p - [-sin(2*x) 1-3*sin(x).^2 sin(2*x)/2])
%
% Example 3:
%  L = chebop(@(x,u,v) [diff(u,2), 0.5*diff(v)+exp(x)]);
%  p = recoverCoeffs(L)
%  norm([p{:}] - [0 0 1 0 0 0 .5])

% Convert to linop if input is a chebop. (But don't overwrite input as it's
% more efficient to evaluate the chebop .op than the linearised .oparray!)
if isa(L, 'chebop')
    L2 = linop(L);
else
    L2 = L;
end

% Initialise:
s = size(L2);                    % Determine the size of the system,
m = L2.diffOrder;                % and the difforder.
x = chebfun('x', L2.domain,2);   % Construct linear function on the domain,
x0 = chebfun(0, L2.domain);      % and the zero function.
p = cell(s);                     % Initialise output.
p0 = L*repmat(x0, 1, s(2));      % Compute non-autonomous component.

% The main routine:
for hh = 1:s(2)                 % Loop over each of the dependent variables.
    x0l = repmat(x0,1,hh-1);    % Set dep vars to the left to zero.
    x0r = repmat(x0,1,s(2)-hh); % Set dep vars to the right to zero.
    p1 = L*[x0l 1+0*x x0r];     % Evaluate all equations for [0 ... 1 ... 0]
    p1 = p1 - p0;               % Subtract non-autonomous compnent.
    for ll = 1:s(1)             % Loop over equations and assign.
        p{ll,hh} = p1{ll};
    end
    xk = x;                            % Update indep var to x.
    for kk = 1:max(m(:,hh))            % Loop over each x^k.
        tmp = L*[x0l xk(:,kk) x0r]-p0; % Evaluate for u = [0 ... x^k ... 0].
        for ll = 1:s(1)                % Loop over each equation.
            if kk > m(ll,hh)           % No coeffs of this order here.
                continue
            end 
            p{ll,hh}(:,kk+1) = tmp{ll}; % Assign the ll-th equation.
            for jj = 1:kk               % Extract the lower-order terms.
                p{ll,hh}(:,kk+1) = p{ll,hh}(:,kk+1) - p{ll,hh}(:,kk+1-jj).*xk(:,jj);
                p{ll,hh}(:,kk+1) = simplify(p{ll,hh}(:,kk+1)); % Simplify.
            end
        end
        xk = [xk x.*xk(:,end)/(kk+1)]; % Update indep var to x^k/k!
    end
end

end

function val = chebValues(k, n, x)
%CHEBBALUES   Return the values of Chebyshev {T0^(k)(x),..Tn^(k)(x)}, x being 
% 1 or -1. 
if ( k == 0 )
    val = x.^((0:n-1).');
else
    [ll, kk] = meshgrid((0:n-1), (0:k-1));
    val = (x).^((1:n).').*prod((ll.^2 - kk.^2)./(2*kk+1), 1).';
end

end
