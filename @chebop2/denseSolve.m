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
            
%            debug:
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

function X = ImposeBoundaryConditions(X, bb, gg, Px, Py, m, n)

% Recombine in the boundary conditions.
cs = size(bb{3},2) + size(bb{4},2);
rs = size(bb{1},2) + size(bb{2},2);

By = [bb{3}.' ; bb{4}.'].'; Gy = [gg{3}.' ; gg{4}.'].';
if ( ~isempty(By) )
    By = Py.' * By; % Gy = Gy * Py;
    X12 = By(1:cs,:).' \ (Gy(rs+1:size(X,2)+rs,:).' - By(cs+1:size(X,1)+cs,:).'*X);
    X = [ X12; X ];
end

Bx = [bb{1}.' ; bb{2}.'].'; Gx = [gg{1}.' ; gg{2}.'].';
if ( ~isempty(Bx) )
    Bx = Px.' * Bx; %Gx = Px.' * Gx;
    X2 = ( Bx(1:rs,:).' \ ( Gx(1:size(X,1),:).' - Bx(rs+1:size(X,2)+rs,:).'*X.') ).' ;
    X = [X2 X];
end

if ( size(X,1) < m )
    X(size(X,1)+1:m, :) = 0;
end

if ( size(X,2) < n )
    X(:, size(X,2)+1:n) = 0;
end

if ( ~isempty(Px) )
    X = X * Px.';
end

if ( ~isempty(Py) )
    X = Py * X;
end

end