function X = denseSolve(N, f, m, n)
% Given a fixed discretisation size, solve the differential equation using
% a dense solver.
%
%  DENSESOLVE(N,F,M,N), returns a solution matrix of values on a M by N
%  Chebyshev grid.

[CC, RHS, bb, gg, Px, Py, xsplit, ysplit] =...
    chebop2.constructDiscretisation(N,f,m,n); % Construct discretisation.

if ( size(CC,1) == 1 )  % rank-1 PDE operator.
    A = CC{1,1}; B = CC{1,2};
    Y = A \ RHS ;
    X = ( B \ Y.').';
    
    X = ImposeBoundaryConditions(X,bb,gg,Px,Py,m,n);
    
elseif ( size(CC,1) == 2 )% rank-2 PDE operator.
    
    A = CC{1,1}; B = CC{1,2}; C = CC{2,1}; D = CC{2,2};
    
    if min(size(bb{1})) > 1 || min(size(bb{2}))>1 ||...
            min(size(bb{3}))>1 || min(size(bb{4}))>1
        xsplit = 0; ysplit = 0;
    end
    
    %xsplit = 0; ysplit = 0; 
    try
        if xsplit || ysplit
            % solve subproblems.
            X = chebop2.BartelsStewart(A,B,C,D,RHS,xsplit,ysplit);
            
%             debug:
%            norm(A * X * B' + C * X * D' - RHS,inf)
        else
            X = lyap(C\A,(B\D).',-(B\(C\RHS).').');
            
        end
    catch
        X = chebop2.BartelsStewart(A,B,C,D,RHS,xsplit,ysplit);% Solve with Sylvester Solver\
    end
    % QZ in BartelsStewart solver can return complex so take real part
    if norm(imag(X),inf) < sqrt(eps)
        X = real(X);
    end
    
    X = ImposeBoundaryConditions(X,bb,gg,Px,Py,m,n);

    if (  any( isnan(X(:)) | isinf(X(:)) )  )
        error('CHEBOP2:MLDIVIDE', 'Solution is not unique.')
    end
    %surf(log10(abs(X)))
else  % rank-k>2 PDE operator.
    % do full n^2 by n^2 matrix kronecker product.
    
    % make massive n^2 by n^2 matrix.
    sz = size(CC{1,1},1)*size(CC{2,1},2);
    if sz > 60^2
        error('CHEBOP2:MLDIVIDE','Solution was unresolved on a 60 by 60 grid.');
    end
    A = spalloc(sz, sz, m*sz + sz);
    for jj = 1:size(CC,1)
        A = A + kron(CC{jj,2},CC{jj,1});
    end
    % vec the rhs.
    b = RHS(:);
    
    X = A \ b;
    if max(size(X)) > 4000
        error('Unresolved for n > 4000.');
    end
    
    X = reshape(X,size(CC{1,1},1),size(CC{2,1},2));
    
    X = ImposeBoundaryConditions(X,bb,gg,Px,Py,m,n);
end

% Debug:
% X = polyval(rot90(X,2));
% x = chebpts(size(X,1));
% y = chebpts(size(X,2));
% [xx,yy]=meshgrid(y,x);
% surf(xx,yy,X,'edgealpha',.5,'facecolor','interp');
%surf(log10(abs(X)))
end