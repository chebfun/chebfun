function varargout = eig(f)
%EIG    Eigenvalues and eigenfunctions of a CHEBFUN2.
%   EIG(F) returns the eigenvalues of F. The number of nonzero eigenvalues
%   returned is at most the rank of the CHEBFUN2. F.domain needs to be
%   square: [a b a b].
%
%   S = EIG(F) returns the eigenvalues of F. S is a vector of 
%   eigenvalues.
%
%   [V,D] = EIG(F) returns the eigenvalues and eigenfunctions of F. V is a quasi-matrices of
%   CHEBFUN objects and D is a diagonal matrix with the eigenvalues
%   on the diagonal.
%
%   If d is an eigenvalue of F with eigenfunction w, then this holds: 
%   \int f(x,y)w(x)dx = w(y). 
%   This can be checked as follows: 
%   [V,D] = eig(f);
%   w = V(:,1); d = D(1,1);
%   fw = chebfun(@(y)sum(chebfun(@(x)f(x,y),w.domain).*w),w.domain);
%   and check that norm(fw-d*w) is O(eps). 
%

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    dom = f.domain; 
if ( (dom(1) ~= dom(3)) || (dom(2) ~= dom(4)) )% domain check
    error('CHEBFUN2:CHEBFUN2:eig:domainerr',...
          'Domain of chebfun2 needs to be in form [a b a b]');
end    

    [u, S, v] = svd( f ); % SVD of chebfun2
    [V, D] = eig( v'*u*S ); % eig(AB) = eig(BA)
    V = u * S * V; % eigenfunctions

    if ( nargout > 1 )
        varargout = { V, D };
    else
        varargout = { diag( D ) };
    end

end
