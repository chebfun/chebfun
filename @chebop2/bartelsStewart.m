function X = bartelsStewart(A, B, C, D, E, xsplit, ysplit)
%BARTELSSTEWART   Solution to generalized Sylvester matrix equation. 
% 
% Computes the solution to the Sylvester equation
%
%         AXB^T + CXD^T = E
%
% by using the Bartels--Stewart algorithm, see 
%
% J. D. Gardiner, A. J. Laub, J. J. Amato, & C. B. Moler, Solution of 
% the Sylvester matrix equation AXB^T+ CXD^T= E, ACM Transactions on 
% Mathematical Software (TOMS), 18(2), 223-231.
% 
% This Bartels--Stewart solver also takes information xsplit, ysplit so
% that if possible it decouples the even and odd modes.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun information.

tol = eps();

% If the RHS is zero then the solution is the zero solution (assuming
% uniqueness).
if ( norm(E) < 10*tol )
    X = zeros(size(E));
    return
end

% If the equation is even/odd in the x-direction then we can split the problem
% into two subproblems. We enforce P and S as upper triangular because they
% are up to rounding errors, and we need to do back substitution with
% them.
if ( ysplit )
    % This is equivalent to qz(full(A),full(C)), but faster.
    [P, S, Q1, Z1] = qzsplit(A, C); 
    P = triu(P); 
    S = triu(S);
else
    [P, S, Q1, Z1] = qz(full(A), full(C));  
    P = triu(P); 
    S = triu(S);
end

% If the PDE is even/odd in the y-direction then we can split (further)
% into double as many subproblems.
if ( xsplit )
    % Faster QZ when even/odd modes decouple in x-direction: 
    [T, R, Q2, Z2] = qzsplit(D, B);
else
    % QZ does not take matrices in sparse format:
    [T, R, Q2, Z2] = qz(full(D), full(B));
end

% Now use the generalised Bartels--Stewart solver found in Gardiner et al.
% (1992).  The Sylvester matrix equation now contains quasi upper-triangular
% matrices and we can do a backwards substitution.

% transform the righthand side.
F = Q1*E*Q2.';

% Solution will be a m by n matrix.
m = size(A, 1); 
n = size(B, 1); 
Y = zeros(m, n);

% Do a backwards substitution type algorithm to construct the solution.
k=n;
PY = zeros(m);
SY = zeros(m);

% Construct columns n,n-1,...,3,2 of the transformed solution.  The first
% column is treated as special at the end.
while k > 1
    % There are two cases, either the subdiagonal contains a zero
    % T(k,k-1)=0 and then it is a backwards substitution, or T(k,k-1)~=0
    % and then we solve a 2x2 system instead.
    
    if ( T(k,k-1) == 0 )
        % Simple case (almost always end up here).
        rhs = F(:,k);
        if ( k < n )
            
            PY(:,k+1) = P*Y(:,k+1);
            SY(:,k+1) = S*Y(:,k+1);
            
            for jj = k+1:n
                rhs = rhs - R(k,jj)*PY(:,jj) - T(k,jj)*SY(:,jj);
            end
            
        end
        
        % find the kth column of the transformed solution.
        Y(:,k) = (R(k,k)*P + T(k,k)*S) \ rhs;
        
        % go to next column
        k = k-1;
        
    else
        % This is a straight copy from the Gardiner et al. paper, and just
        % solves for two columns at once. (works because of
        % quasi-triangular matrices.
        
        % Operator reduction.
        rhs1 = F(:,k-1);
        rhs2 = F(:,k);
        
        for jj = k+1:n
            yj = Y(:,jj);
            rhs1 = rhs1 - R(k-1,jj)*P*yj - T(k-1,jj)*S*yj;
            rhs2 = rhs2 - R(k,jj)*P*yj - T(k,jj)*S*yj;
        end
        
        % 2 by 2 system.
        SM = zeros(2*n);
        up = 1:n;
        down = n+1:2*n;
        
        SM(up,up) = R(k-1,k-1)*P + T(k-1,k-1)*S;
        SM(up,down) = R(k-1,k)*P + T(k-1,k)*S;
        SM(down,up) = R(k,k-1)*P + T(k,k-1)*S;
        SM(down,down) = R(k,k)*P + T(k,k)*S;
        
        % Permute the columns and rows: 
        Spermuted = zeros(2*n);
        Spermuted(1:2:2*n,1:2:2*n) = SM(1:n,1:n); 
        Spermuted(2:2:2*n,2:2:2*n) = SM(n+1:2*n,n+1:2*n); 

        % Solve.
        UM = Spermuted \ [rhs1; rhs2];
        
        Y(:,k-1) = UM(up); 
        Y(:,k) = UM(down);

        PY(:,k) = P*Y(:,k);
        PY(:,k-1) = P*Y(:,k-1);
        SY(:,k) = S*Y(:,k); 
        SY(:,k-1) = S*Y(:,k-1);
        
        % We solved for two columns so go two columns further.
        k=k-2;
        
    end
    
end

if ( k == 1 )
    % Now we have just the first column to compute.
    rhs = F(:,1);
    PY(:,2) = P*Y(:,2);
    SY(:,2) = S*Y(:,2);
    for jj = 2:n
        rhs = rhs - R(1,jj)*PY(:,jj) - T(1,jj)*SY(:,jj);
    end
    Y(:,1) = (R(1,1)*P + T(1,1)*S) \ rhs;
end

% We have now computed the transformed solution so we just transform it
% back.
X = Z1*Y*Z2.';

end

function [P, S, Q1, Z1] = qzsplit(A, C)
%QZSPLIT   A faster qz factorisation for problems that decouple.
%
% This is equivalent to standard qz, except we take account of symmetry to
% reduce the computational requirements of the QZ factorisation.

% Do the QZ by splitting the problem into two subproblems. 
A = full(A); 

A1 = A(1:2:end,1:2:end); 
C = full(C); 
C1 = C(1:2:end,1:2:end);
[P1, S1, Q1, Z1] = qz(A1, C1);

A2 = A(2:2:end,2:2:end); 
C2 = C(2:2:end,2:2:end);
[P2, S2, Q2, Z2] = qz(A2, C2);

[P, S, Q1, Z1] = reform(P1, P2, S1, S2, Q1, Q2, Z1, Z2);

end

function [P, S, Q, Z] = reform(P1, P2, S1, S2, Q1, Q2, Z1, Z2)
%REFORM   Recombine subproblems to form the QZ factorization. 

% Initialise all the variables. 
hf1 = size(P1, 1);
n = 2*hf1 - 1;
P = zeros(n);
S = zeros(n);
Q = zeros(n);
Z = zeros(n);

% Push the subproblem back together.
P(1:hf1,1:hf1) = P1; 
P(hf1+1:end,hf1+1:end) = P2;

S(1:hf1,1:hf1) = S1; 
S(hf1+1:end,hf1+1:end) = S2;

Q(1:hf1,1:2:end) = Q1; 
Q(hf1+1:end,2:2:end) = Q2;

Z(1:2:end,1:hf1) = Z1; 
Z(2:2:end,hf1+1:end) = Z2;

end
