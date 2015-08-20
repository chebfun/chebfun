function [CC, rhs, bb, gg, Px, Py, xsplit, ysplit] = discretize(N, f, m, n, flag)
%DISCRETIZE   Given a CHEBOP2, this function converts the problem to one of the form
%
%  sum_i  kron(A_i,B_i)
%
% and computes the discretisation as a cell array: 
%
%  {{
%
%       A_1 , B_1
%       A_2 , B_2
%        .  ,  .
%        .  ,  .
%        .  ,  .
%       A_k , B_k
%
%                   }}
%
%    INPUTS: 
%      N = PDE (CHEBOP2). 
%      f = forcing term (CHEBFUN2).
%      m = discretization size in 1st variable.
%      n = discretization size in 2nd variable.
%      flag = 0 (default) means assigned boundary conditions, flag = 1
%      means do not assign boundary conditions. 
%
%    OUTPUTS:
%      CC = cell array of matrices storing the terms in the matrix
%      equation.
%      rhs = matrix discretizing forcing term. 
%      bb = cell array storing the discretized linear constraints. 
%      gg = cell array storing the discretized nonhomogeneous part of the 
%      constraints. 
%      Px, Py = store permutation matrix to ensure bcs are linear
%      dependent.
%      XSPLIT, YSPLIT = 0 if subproblems cannot be formed. XSPLIT = 1 if
%      even and odd modes decouple in 1st variable. YSPLIT = 1 if even and
%      odd modes decouple in 2nd variable. 
%
% Returns RHS with degrees of freedom removed and bb which stores
% elminated boundary conditions, gg eliminated boundary rows.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs.
if nargin < 5
    flag = 0;
end

%%
% Get information of PDE and prefs.
A = N.coeffs;
rect = N.domain;
prefs = chebfunpref(); 
tol = prefs.techPrefs.eps; 
xorder = N.xorder;
yorder = N.yorder;

%%
% Check if the PDO was given as a variable coefficient PDO with the notation
% @(x,y,u), but is actually a constant coefficient PDO. 
if ( iscell( A ) ) 
    doesNotDependOnXorY = all(all(cellfun(@isnumeric, N.coeffs))); 
    if ( doesNotDependOnXorY ) 
        A = cell2mat( A ).'; 
    end
end

%%
% Convert matrix of coefficients to a discretization for the PDE using the
% singular value decomposition.  We find the rank of the PDE operator and
% then use the optimal low rank expansion of the operator as a way to
% discretise the PDE.
if (  isempty(N.V) || isempty(N.U) )
    if ( iscell(A) )
        counter = 1;
        U = cell(size(A, 1), 1); 
        V = cell(size(A, 2), 1);
        for jj = 1:size(A, 1)
            for kk = 1:size(A, 2)
                a = A{jj,kk};
                if ( isa(a, 'chebfun2') )
                    if ( abs(vscale(a)) > tol )
                        [C, D, R] = cdr(a);
                        for col = 1:size(C, 2)
                            U{jj,counter} = C(:,col)*D(col,col);
                            V{kk,counter} = R(:,col);
                            counter = counter + 1;
                        end
                    end
                elseif ( isa(a,'double') && ( isempty(a) || abs(a) > tol ) )
                    U{jj,counter} = a; 
                    V{kk,counter} = 1;
                    counter = counter + 1;
                end
            end
        end
        rk = size(U, 2);
        S = diag(ones(rk, 1));
        na = size(U, 1);
        nb = size(V, 1);
    else
        % Compute the SVD of the coefficient matrix.
        [U, S, V] = svd(A.');
        % Find the rank of A, which is also the rank of the PDE operator and
        % construct the low rank expansion for A.
        rk = find(diag(S) > tol, 1, 'last');
        U = U(:,1:rk);
        S = S(1:rk,1:rk);
        V = V(:,1:rk);
        [na, nb] = size(A.');
    end
else
    rk = size(N.S, 2);
    U = N.U;
    V = N.V;
    S = N.S;
    na = size(U, 1); 
    nb = size(V, 1);
end
% LEFT = zeros(m); RIGHT = zeros(n);

% Construct the discretisation in matrix equation form.
CC = cell(rk,2);

for jj = 1 : rk
    
    RIGHT = unconstrainedMatrixEquation(V, jj, n, xorder, rect(1:2)); % jjth term on the right.
    LEFT = unconstrainedMatrixEquation(U, jj, m, yorder, rect(3:4)); % jjth term on the left.
    
    % Balance out the scaling from the singular value. This does slightly
    % improve the accuracy.
    singvalue = sqrt(S(jj,jj));
    CC{jj,2} = singvalue * RIGHT;
    CC{jj,1} = singvalue * LEFT;
    
end

%%
% Test to see if we can solve subproblems. This checks if the PDE operator
% contains differential terms of the same parity.
ysplit = 0; xsplit=0;
if ~iscell(U) && ~iscell(V)
    emask = 1:2:na;
    omask = 2:2:na;
    if ( min( norm(U(emask,:)), norm(U(omask,:)) ) < 10*tol )
        ysplit = 1;
    end
    emask = 1:2:nb;
    omask = 2:2:nb;
    if ( min( norm(V(emask,:)), norm(V(omask,:)) ) < 10*tol )
        xsplit = 1;
    end
end

%%
% We have a discretisation for the PDE operator, now let's find a
% discretisation for the boundary conditions.

% If no boundary conditions is prescribed then make it empty.
bcLeft = []; 
leftVal = [];
bcRight = []; 
rightVal = [];
bcUp = []; 
upVal = [];
bcDown = []; 
downVal = [];

if ( ~isempty(N.lbc) ) % Left boundary conditions.
    [bcLeft, leftVal] = chebop2.constructBC(N.lbc, -1, m, n, rect(3:4), rect(1:2), xorder);
end
if ( ~isempty(N.rbc) ) % Right boundary conditions.
    [bcRight, rightVal] = chebop2.constructBC(N.rbc, 1, m, n, rect(3:4), rect(1:2), xorder);
end
if ( ~isempty(N.ubc) ) % Top boundary conditions.
    [bcUp, upVal] = chebop2.constructBC(N.ubc, 1, n, m, rect(1:2), rect(3:4), yorder);
end
if ( ~isempty(N.dbc) ) % Bottom boundary conditions.
    [bcDown, downVal] = chebop2.constructBC(N.dbc, -1, n, m, rect(1:2), rect(3:4), yorder);
end

%%

% For the down and up BCs we have B^TX = g^T.
By = [ bcUp.'; bcDown.' ];
Gy = [ upVal.'; downVal.' ];
[By, Gy, Py] = canonicalBC(By, Gy);

% For the left and right BCs we have X*B = g. We do the LU to B^T.
Bx = [ bcLeft.'; bcRight.' ];
Gx = [ leftVal.'; rightVal.' ];
[Bx, Gx, Px] = canonicalBC(Bx, Gx);
Bx = Bx.';
Gx = Gx.'; % Now transpose so that X*B = g;

%% 
% Construct the RHS of the Sylvester matrix equation.

% Complete the RHS (part of the RHS could have been in the operator): 
if ( ~isempty(N.rhs) && ( isa(N.rhs, 'chebfun2') || isa(N.rhs, 'double') ) )
    f = f + N.rhs; 
end
E = zeros(m, n);
[n2, n1] = length(f);
F = chebcoeffs2(f);

% Map the RHS to the right ultraspherical space.
lmap = ultraS.convertmat(n1, 0, yorder-1);
rmap = ultraS.convertmat(n2, 0, xorder-1);
F = lmap * F * rmap.';

% Place those coefficients of the forcing function onto the RHS.
n1 = min(n1, m); 
n2 = min(n2, n); 
E(1:n1,1:n2) = F(1:n1,1:n2);

if ( ~flag ) % Impose boundary conditions.
    
    % Use the eliminated boundary condition to place zeros in the columns of
    % the matrix equation discretization. There are rk columns to zero out.
    
    for jj = 1:rk  % For term in the matrix equation.
        [C, E] = zeroDOF(CC{jj,1}, CC{jj,2}, E, By, Gy);
        CC{jj,1} = C;
        [C, E] = zeroDOF(CC{jj,2}, CC{jj,1}, E.', Bx.', Gx.');
        CC{jj,2} = C; 
        E = E.';
    end
    
    % Remove degrees of freedom.
    nn = n - max(xorder, yorder);
    mm = m - max(xorder, yorder);
    df1 = max(0, xorder - yorder);
    df2 = max(0, yorder - xorder);
    for jj = 1:rk
        CC{jj,1} = CC{jj,1}(1:mm, yorder+1:m-df1);
        CC{jj,2} = CC{jj,2}(1:nn, xorder+1:n-df2);
    end
    % Truncation of righthand side.
    rhs = E(1:mm, 1:nn);
    
else
    rhs = E;
end

% Pass back the eliminated boundary conditions.
bb = {bcLeft bcRight bcUp bcDown};
gg = {leftVal rightVal upVal downVal};

%% 
% Check boundary continunity conditions.

% Check BCs at corners:
allbc = 0;
if ( ~isempty(bcUp) && ~isempty(upVal) && ~isempty(bcRight) && ~isempty(rightVal) )
    if ( norm(rightVal(end-4:end),inf) < sqrt(tol) && norm(upVal(end-4:end),inf) < sqrt(tol) )
        allbc = allbc + norm(upVal.'*bcRight - bcUp.'*rightVal);
    end
end
if ( ~isempty(bcUp) && ~isempty(upVal) && ~isempty(bcLeft) && ~isempty(leftVal) )
    if ( norm(leftVal(end-4:end),inf) < sqrt(tol) && norm(upVal(end-4:end),inf) < sqrt(tol) )
        allbc = allbc + norm(upVal.'*bcLeft - bcUp.'*leftVal);
    end
end
if ( ~isempty(bcDown) && ~isempty(downVal) && ~isempty(bcRight) && ~isempty(rightVal) )
    if ( norm(rightVal(end-4:end),inf) < sqrt(tol) && norm(downVal(end-4:end),inf) < sqrt(tol) )
        allbc = allbc + norm(downVal.'*bcRight - bcDown.'*rightVal);
    end
end
if ( ~isempty(bcDown) && ~isempty(downVal) && ~isempty(bcLeft) && ~isempty(leftVal) )
    if ( norm(leftVal(end-4:end),inf)<sqrt(tol) && norm(downVal(end-4:end),inf)<sqrt(tol) )
        allbc = allbc + norm(downVal.'*bcLeft - bcDown.'*leftVal);
    end
end

% [TODO]: Should there be an error if the compatibility conditions do not 
% match?
%
% if allbc >= 100*sqrt(tol)
%     s = sprintf('Boundary conditions differ by %1.4f', allbc');
%     warning('CHEBFUN:CHEBOP2:discretize:BCs', s)
% end

end

function B = unconstrainedMatrixEquation(ODE, jj, n, order, dom)
% Construct the unconstrained Matix Equation. Adding in the constraints later.

B = spalloc(n, n, 3*n);
for kk = 1:size(ODE, 1)
    
    % Get conversion and differentiation matrices: 
    S = ultraS.convertmat(n, kk-1, order-1);
    D = ((2./diff(dom))^(kk-1)) * ultraS.diffmat(n, kk-1);
    
    if ( iscell(ODE(kk,jj)) && isa(ODE{kk,jj}, 'chebfun') )
        % Variable coefficient term: 
        c = ODE{kk,jj}.coeffs;        
        M = ultraS.multmat(n, c, kk-1); 
        A = S * M * D;
        
    elseif ( iscell(ODE(kk,jj)) && ~isempty(ODE{kk,jj}) )
        % Constant coefficient term in a variable coefficient ODE:
        A = ODE{kk,jj}.* S * D;
        
    elseif ( isa(ODE(kk,jj),'double') )
        % Constant coefficient term in a constant coefficient ODE:
        A = ODE(kk,jj).* S * D;
        
    else
        % Empty cell in array, so no term: 
        A = zeros(n);
        
    end
    
    % Form the ODE operator: 
    B = B + A;
    
end

end

function [B, G, P] = canonicalBC(B, G)
%CANONICALBC   Form a linear combintation of the boundary conditions 
%so that they can be used for imposing on the PDE. 

P = nonsingularPermute(B);
B = B*P;
[L, B] = lu(B); 
G = L \ G;

% Scale so that B is unit upper triangular.
if ( min(size(B)) > 1 )
    D = diag(1./diag(B));
elseif ( ~isempty(B) )
    D = 1./B(1,1);
else
    D = []; % No boundary conditions.
end
B = D*B; 
G = D*G;

end

function P = nonsingularPermute(B)
%NONSINGULARPERMUTE   Permute the columns of B to ensure that the principal
%m*m submatrix of B is nonsingular, where m = size(B, 1).
%
% Note: This is needed for solving the matrix equations with linear
% constraints, see DPhil thesis of Alex Townsend (section 6.5).

m = size(B, 1);
k = 1;

% [TODO]: improve this check.
% Try each mxm block in a linear fashion: 
while ( rank(B(:,k:m+k-1)) < m )
    k = k+1;
    if ( m+k > size(B, 2) )
        error('CHEBFUN:CHEBOP2:discretize:nonsingularPermute:BCs', ...
            'Boundary conditions are linearly dependent.');
    end
end

P = speye(size(B, 2));
P = P(:,[k:m+k-1, 1:k-1, m+k:end]);

end

function [C1, E] = zeroDOF(C1, C2, E, B, G)
%ZERODOF   Eliminate so degrees of freedom in the matrix equation can be
%removed.

for ii = 1:size(B, 1) % For each boundary condition, zero a column.
    for kk = 1:size(C1, 1)
        if ( abs(C1(kk,ii)) > 10*eps )
            c = C1(kk, ii); % Constant required to zero entry out.
            C1(kk,:) = C1(kk,:) - c*B(ii,:);
            E(kk,:) = E(kk,:) - c*G(ii,:)*C2.';
        end
    end
end

end
