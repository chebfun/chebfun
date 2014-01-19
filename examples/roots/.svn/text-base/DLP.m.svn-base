function [X, Y, M] = DLP(AA, v, a, b, c)
pwd
%DLP constructs the DL pencil with ansatz vector v. 
% 
% [X,Y] = DLP(AA,V,A,B,C) returns the DL pencil with the orthogonal
%    basis defined by the recurrence relations A,B,C.

[n,m] = size(AA); k=m/n-1; s=n*k;              % matrix size and degree  
%M = spdiags([a b c],[0 n 2*n],s,s+n);          % multiplication matrix
M = spdiags([a b c],[0 1 2],k,k+1); M = kron(M,eye(n)); % multiplication matrix

S = kron(v,AA);
for j=0:k-1, jj=n*j+1:n*j+n; AA(:,jj) = AA(:,jj)';end % block transpose
T = kron(v',AA'); R = M'*S-T*M;                % construct RHS

% Bartel-Stewart algorithm on M'Y+YM=R, M is upper triangular. 
X = zeros(s); Y=X; ii=n+1:s+n; nn=1:n;         % useful indices
Y(nn,:)=R(nn,ii)/M(1); X(nn,:)=T(nn,:)/M(1);   % first column of X and Y
Y(nn+n,:)=(R(nn+n,ii)-M(1,n+1)*Y(nn,:)+Y(nn,:)*M(:,n+1:s+n))/M(n+1,n+1);
X(nn+n,:)=(T(nn+n,:)-Y(nn,:)-M(1,n+1)*X(nn,:))/M(n+1,n+1); % 2nd columns

for i = 3:k                                    % backwards substitution
    ni=n*i; jj=ni-n+1:ni; j0=jj-2*n; j1=jj-n;  % useful indices
    M0=M(ni-2*n,ni); M1=M(ni-n,ni); m=M(ni,ni);% consts of 3-term recurr
    Y0=Y(j0,:); Y1=Y(j1,:); X0=X(j0,:); X1=X(j1,:);% variables in 3-term
    Y(jj,:)=(R(jj,ii)-M1*Y1-M0*Y0+Y1*M(:,n+1:s+n))/m; 
    X(jj,:)=(T(jj,:)-Y1-M1*X1-M0*X0)/m;        % use Y to solve for X
end