function  [ZZ, DD, YY] = fiadi(A,B,U,S,V,varargin)
%
% Applies factored-independent ADI to solve AX - XB = U*S*V', where F =U*S*V' 
% is an (approximate) SVD of F and F has rapidly decaying singular values.
%
% The equation AX-XB = F is split into N equations, 
% 
%   AX_i + X_iB = F_i, 
% 
% where X_1 + ... + X_N = X, F_1 + ... + F_N = F. Factored ADI is then 
% independently applied to each equation. 
% 
% fiadi(A, B, U, S, V, P, Q,tol) applies FIADI with the nonzero entries in
% P(i,:), Q(i, :) used as the set of shift parameters applied to the ith 
% equation above. (See getshifts_fiadi, getshifts_fismith). The solution is 
% compressed using the relative tolerance tol. If tol is not specified, 
% it is automatically set to machine precision.
%
% fiadi(A, B, U, S, V, P, Q, 'nocompress') is the same as above, except that the
% uncompressed solution is returned.
%
% fiadi(A, B, U, S, V) If A and B have eigenvalues contained in real,
% disjoint intervals, then the fiadi shift parameters are computed
% automatically and the solution is computed to approximately machine precision.
%
% fiadi(A, B, U, S, V, tol) Same as above except that the solution is 
% computed to a relative accuracy specified by tol. 
%
% References:
% [1] Townsend, Alex, and Heather Wilber. 
% "On the singular values of matrices with high displacement rank." 
% arXiv preprint arXiv:1712.05864 (2017).

% code written by Heather Wilber (heatherw3521@gmail.com)

%%
%
[m, ~]= size(U);
[n,~] = size(V); 
compute_shifts =0; 

s = diag(S); 
if ~isreal(s)
    s = real(s); 
    warning('ADI:FIADI',...
    'Nonzero imaginary part of approximate singular values detected. Setting singular values to real part only.')
end
S = diag(s); 
U = U*sqrt(S); 
V = V*sqrt(S); 
rnk = length(s); 

% parse inputs: 

if isempty(varargin) 
    comp = 1;
    tol = 1e-16; 
    compute_shifts = 1; 
elseif numel(varargin)==1 %no shift params, either specified tol or nocompress
    compute_shifts = 1; 
    if strcmp(varargin{1}, 'nocompress')
        comp = 0; 
        tol = eps;
    else
        comp = 1; 
        tol = varargin{1}; 
    end
elseif numel(varargin)==2 %shift parameters given but no compression tolerance param.
        P = varargin{1}; 
        Q = varargin{2}; 
        comp = 1; 
        tol = eps;
elseif numel(varargin)==3 %shift parameters and tolerance for compression.
        P = varargin{1}; 
        Q = varargin{2}; 
        if strcmp(varargin{3}, 'nocompress') 
        comp = 0; 
        else 
            tol = varargin{3}; 
            comp = 1; 
        end
else %not sure what to do
    error('adi:fiadi: input invalid.')
end
    

%user wants shift parameters computed
if compute_shifts == 1 
    %find intervals where eigenvalues live: 
    a = eigs(A, 1, 'SM'); 
    b = eigs(A, 1, 'LM'); 
    c = eigs(B, 1, 'SM'); 
    d = eigs(B, 1, 'LM'); 
    % determine if the eigenvalues have a complex part. if so, abandon. 
    if any(abs(imag([a b c d])) > 1e-10)
        error('ADI:fadi:cannot automatically compute shift parameters unless the eigenvalues of A and B in AX - XB = F are contained in real, disjoint intervals')
    end
    % check if intervals overlap
    I1 = [min(a,b) max(a,b)]; 
    I2 = [min(c,d) max(c,d)]; 
    if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
        error('ADI:fadi:cannot automatically compute shift parameters unless the eigenvalues of A and B in AX - XB = F are contained in real, disjoint intervals.')
    end
    [P,Q] = getshifts_fiadi([I1 I2],A, B, U, S, V, tol); 
end

%%
% 
% we set up a block parameter associated with how each F_i is formed. 
% The specific strategy we use here is very simple: each F_i = U(:,i)*S(i,i)*V(:,i)', 
% and then we group together RHSs that use the same number of shift parameters. 

%set up blocks = [nrhs, nitr]; 
% when the blocks are too big, it can make compression expensive. So we set a 
% max block size parameter. 

count = sum(abs(P)>0,2);
mxsz = floor(.3*max(n,m)); 

nr = 0;
nc = 0; 
ncold = 0; 
blocks = []; 
ncidx = [];
while(nc < rnk ) 
ni = count(nc+1);
nc = find(count==ni, 1,'last');
nr = min( (nc-ncold), mxsz);
nc = ncold + nr; 
blocks = [blocks; nr ni];
ncold = nc;
ncidx = [ncidx; nc]; 
end

%eliminate last block if zero ADI steps: 
zidx = find(blocks(:,2)==0, 1, 'first');
if ~isempty(zidx)
blocks = blocks(1:zidx-1, :); 
ncidx = ncidx(1:zidx-1); 
end

%adjust P and Q: 
P = P(ncidx,:); 
Q = Q(ncidx, :); 
Nblk = length(blocks(:,1)); %total number of blocks to compute.


%% perform fadi on each block
r = 0; 
rold=0;
YY =[];
ZZ =[];
DD = []; 
for k = 1:Nblk
    rold = r+rold; %pick out the correct RHS block
    r = blocks(k,1);
    p = P(k,:); 
    q = Q(k,:); 
    p = p(~(p==0)); 
    q = q(~(q==0));
    C1 = U(:,rold+1:r+rold);
    R1 = V(:,rold+1:r+rold);
      
    % Do fadi on the block
    [Z, D, Y] = fadi(A, B, C1,R1,p,q);
         
    D = diag(D); 
    if comp == 1  
    % Compress after each block
        ZZ = [ZZ Z]; 
        YY = [YY Y];        
        DD = [DD; D];
       [ZZ, DD, YY] = compression(ZZ, DD, YY, tol, 1);
       DD = diag(DD); 
    else 
        ZZ = [ZZ Z]; 
        YY = [YY Y];        
        DD = [DD; D];   
    end  
end
DD = diag(DD);
end





