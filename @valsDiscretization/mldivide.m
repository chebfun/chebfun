function [v, disc] = mldivide(disc, A, b)
%MLDIVIDE  Overloads the default A\b for a discrete linear algebra problem. 
%   V = MLDIVIDE(DISC, A, B) solves V = A\B If a valid factorization of A
%   is stored in DISC, it is used. 
%
%   [V, DISC] = MLDIVIDE(DISC, A, B) stores the computed LU factors in DISC.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isFactored(disc) )
    
    % Use the existing factorization. 
    [L, U, p, s] = deal(disc.mldivideData{:});
    
else
    
    % If there are no usable factors, find them:
    s = 1./ max(1, max(abs(A), [], 2) );  % Row scaling to improve accuracy
    A = bsxfun(@times, s, A);
    [L, U, p] = lu(A, 'vector');
    
    % Store factors:
    disc.mldivideData = {L, U, p, s};
    
end

% Solve the system:
sb = s.*b;
v = U \ ( L \ sb(p) );

end
