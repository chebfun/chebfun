function diffMat = spdiffmat(n,k,varargin)
%DIFFMAT(N,K), computes the kth order US derivative matrix 

if(n==0), diffMat=[]; return; end
if(k==0), diffMat=speye(n); return; end

diffMat = spdiags((0:n-1)',1,n,n);
for s = 1:k-1
    diffMat = spdiags(2*s*ones(n,1),1,n,n)*diffMat;
end

if nargin == 3
    dom = varargin{1}; 
    diffMat = ((2./diff(dom))^k)*diffMat; 
end

end