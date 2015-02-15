function varargout = chol( f, varargin )
%CHOL    Cholesky factorization of a chebfun2. 
%
% R = CHOL( F ), if F is a nonnegative definite chebfun2 then this 
% returns an upper triangular quasimatrix so that R'*R is a
% decomposition of F. If F is not nonnegative definite then an error is thrown.
%
% L = CHOL(F, 'lower'), if F is a nonnegative definite chebfun2 then this
% produces a lower triangular quasimatrix so that L*L' is a decomposition of F.
% If F is not nonnegative definite then an error is thrown. 
% 
% [R, p] = CHOL( F ), with two outputs never throwns an error message. If F is
% nonnegative definite then p is 0 and R is the same as above. If F is 
% symmetric but negative definite or semidefinite then p is a positive 
% integer such that R has p columns and R'*R is a rank p nonnegative definite chebfun2 that approximates F. 
% This is particular useful when F is nonnegative definite, but rounding errors
% have perturbed it to be semidefinite. 
%
% [L, p] = CHOL(F, 'lower') same as above but the first argument is lower 
% triangular. 
%
% For more information about the factorization: 
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, submitted, 2014. 
%
% See also LU, and QR. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    varargout = cell(1, nargout);
    return
end

% Is the chebfun2 on a square domain?: 
if ( ~domainCheck(f.cols, f.rows) )
    error('CHEBFUN:CHEBFUN2:chol:domain', ...
        'Chebfun2 is not on a square domain.');
end

% Get rank of f: 
k = length( f ); 

% All the pivots should be on the y = x line. 
PivLoc = f.pivotLocations; 

if ( isempty( PivLoc ) ) 
    % For some reason (probably because f was made with sampling data, rather
    % than a handle) f didn't have any pivot information. Make it now: 
    f = chebfun2( @(x, y) feval(f,x,y), f.domain);
    varargout = {chol( f )}; 
    return
end

% Find the first pivot location is off-diagonal: 
Diagk = find( PivLoc(:,1) ~= PivLoc(:,2), 1, 'first'); 
k = min( k, Diagk ); 

% Check that the pivots are positive:  
Posk = find( pivots( f ) < 0, 1, 'first');
k = min(k, Posk); 

% Were all the pivots nonnegative and the locations on the diagonal?
posdef = ( k == length( f ) );


% Return an error if the function is not nonnegative definite: 
if ( nargout < 2 && ~isempty( posdef ) && ~posdef )
    error('CHEBFUN:CHEBFUN2:chol:definite', ...
        'Chebfun2 is not nonnegative definite.');
end

% Get the CDR decomposition (already computed by constructor):
[C, D, R] = cdr( f );

% Return an error if the function is not symmetric: 
dom = f.domain; 
r = 0.0192475;   % arbitrary point in [-1,1] 
s = -.34756987;  % arbitrary point in [-1,1]
r = diff(dom(1:2))*r - dom(1); 
s = diff(dom(3:4))*s - dom(3); 
symTest = ( abs(feval(f,r,s) - feval(f,s,r)) < 10*max(epslevel(f.cols),epslevel(f.rows)) );
if ( ~symTest )
    error('CHEBFUN:CHEBFUN2:chol:symmetric', ...
        'The chebun2 must be a symmetric function.');
end

% How many terms are posdef: 
p = k; 
if ( isempty( p ) )
    p = length( f ); 
end
% Extract out posdef part: 
C = C( :, 1:p ); 
D = D( 1:p, 1:p ); 
R = R( :, 1:p ); 

% Balance scaling:
C = C * sqrt( D ); 
R = R * sqrt( D );

% Output to user: 
if ( ( nargout < 2 ) && ( nargin == 1 ) )
    varargout = { R.' };
elseif ( ( nargout < 2 ) && strcmpi( varargin{ 1 }, 'lower') )
    varargout = { C }; 
elseif ( ( nargout == 2 ) && ( nargin == 1 ) )
    varargout = { R.', p }; 
elseif ( ( nargout == 2 ) && strcmpi( varargin{ 1 }, 'lower') )
    varargout = { C, p };
end

end
