function X = denseSolve(N, f, m, n)
%DENSESOLVE   Given a fixed discretisation size, solve the differential equation using
%a dense solver.
%   X = DENSESOLVE(N, F, M, N), returns a solution matrix X of values on a M by N
%  Chebyshev grid.
%
% For further details about the PDE solver, see: 
% A. Townsend and S. Olver, The automatic solution of partial differential
% equations using a global spectral method, in preparation, 2014.
% 
% Warning: This PDE solver is an experimental new feature. It has not been
% publicly advertised. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Construct discretisation for PDE:
[CC, RHS, bb, gg, Px, Py, xsplit, ysplit] =...
    chebop2.discretize(N, f, m, n); 

% Rank-1 PDE operator.
if ( size(CC, 1) == 1 )  
    % We can do this one in our sleep: 
    A = CC{1,1}; 
    B = CC{1,2};
    Y = A \ RHS;
    X = (B \ Y.').';
    X = imposeBoundaryConditions(X, bb, gg, Px, Py, m, n);
    
% Rank-2 PDE operator.
elseif ( size(CC, 1) == 2 )
    % Extract out generatlized matrix equation: 
    A = CC{1,1}; 
    B = CC{1,2}; 
    C = CC{2,1}; 
    D = CC{2,2};
    
    % Don't solve for subproblems if we have lots of bcs on one edge:
    if ( min(size(bb{1})) > 1 || min(size(bb{2}))>1 ||...
            min(size(bb{3}))>1 || min(size(bb{4}))>1 )
        xsplit = 0; 
        ysplit = 0;
    end
    
    % xsplit and ysplit tell the solver if it is possible to solve
    % subproblems. Wrap in a try-catch statement in case LYAP is not on the
    % MATLAB path. 
    try
        if ( xsplit || ysplit )
            % Solve AXB^T + CXD^T = RHS, by doing subproblems: 
            X = chebop2.bartelsStewart(A, B, C, D, RHS, xsplit, ysplit);       
        else
            % Solve AXB^T + CXD^T = RHS with MATLAB LYAP command:
            X = lyap(C\A,(B\D).', -(B\(C\RHS).').');
        end
    catch
        % Solve AXB^T + CXD^T = RHS:
        X = chebop2.bartelsStewart(A, B, C, D, RHS, xsplit, ysplit);
    end
    
    % QZ in bartelsStewart solver can return small complex part. Make it 
    % real.
    if ( norm(imag(X),inf) < sqrt(eps) )
        X = real(X);
    end
    
    % Impose the boundary conditions: 
    X = imposeBoundaryConditions(X, bb, gg, Px, Py, m, n);

    % Quick check that nothing was singular: 
    if ( any( isnan(X(:)) | isinf(X(:)) ) )
        error('CHEBFUN:CHEBOP2:denseSolve:nanInf', 'Nonunique solution to PDE.')
    end
    
% Rank-k > 2 PDE operator.
else  
    % Do full n^2 by n^2 matrix kronecker product.
    % Make massive mn by mn matrix.
    sz = size(CC{1,1}, 1) * size(CC{2,1}, 2);
    if sz > 65^2
        error('CHEBFUN:CHEBOP2:denseSolve:unresolved1', ...
            'Solution was unresolved on a 60 by 60 grid.');
    end
    
    % Form the massive mn by mn matrix.
    A = spalloc(sz, sz, m*sz + sz);
    for jj = 1:size(CC, 1)
        A = A + kron(CC{jj,2}, CC{jj,1});   
    end
    
    % Vectorize the rhs.
    b = RHS(:);
    
    % Solve linear system: 
    X = A \ b;
    
    % Avoid filling the memory: 
    if ( max(size(X)) > 4000 )
        error('CHEBFUN:CHEBOP2:denseSolve:unresolved2', ...
            'Unresolved for n > 4000.');
    end
    
    % Unvectorize.
    X = reshape(X, size(CC{1,1}, 1), size(CC{2,1}, 2));
    
    % Impose linear constraints:
    X = imposeBoundaryConditions(X,bb,gg,Px,Py,m,n);
    
end

end

function X = imposeBoundaryConditions(X, bb, gg, Px, Py, m, n)
% This command imposes the boundary condition on the solution. 
% 
% X = solution with conditions, 
% bb = cell array of linear constraints,
% gg = cell array of data, 
% Px = permutation matrix, identifying the nonsingular block in left/right bcs, 
% Py = permutation matrix, identifying the nonsingular block in top/bottom bcs, 
% m = discretization size in 2nd variable,
% n = discretization size in 1st variable.

% Recombine in the boundary conditions.
cs = size(bb{3}, 2) + size(bb{4}, 2);
rs = size(bb{1}, 2) + size(bb{2}, 2);
By = [ bb{3}.'; bb{4}.' ].'; 
Gy = [ gg{3}.'; gg{4}.' ].';

if ( ~isempty(By) )
    By = Py.' * By; % Gy = Gy * Py;
    X12 = By(1:cs,:).' \ (Gy(rs+1:size(X, 2)+rs,:).' - By(cs+1:size(X, 1)+cs,:).'*X);
    X = [ X12; X ];
end

Bx = [bb{1}.' ; bb{2}.'].'; Gx = [gg{1}.' ; gg{2}.'].';

if ( ~isempty(Bx) )
    Bx = Px.' * Bx; %Gx = Px.' * Gx;
    X2 = (Bx(1:rs,:).' \ ( Gx(1:size(X, 1),:).' - Bx(rs+1:size(X, 2)+rs,:).'*X.')).';
    X = [ X2, X ];
end

% Pad with zeros coefficients:
if ( size(X, 1) < m )
    X(size(X, 1)+1:m,:) = 0;
end
if ( size(X, 2) < n )
    X(:,size(X, 2)+1:n) = 0;
end

% Permute back to original solution:
if ( ~isempty(Px) )
    X = X * Px.';
end
if ( ~isempty(Py) )
    X = Py * X;
end

end
