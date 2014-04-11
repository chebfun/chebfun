function [v, disc] = mldivide(disc, A, b)
%MLDIVIDE  Overloads the default A\b for a discrete linear algebra problem. 
%   If a valid factorization of A is stored in DISC, it is used. Otherwise the
%   factors are found and stored in the DISC output as well. The first output is
%   always v = A\b. 

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( isFactored(disc) )
    
    % Use the existing factorization. 
    [L, U, s] = deal(disc.mldivideData{:});
    
else
    
    % If there are no usable factors, find them:
    s = 1./ max(1, max(abs(A), [], 2) );  % row scaling to improve accuracy
    A = bsxfun(@times, s, A);
    [L, U] = lu(A);
    disc.mldivideData = {L, U, s};
    
end

% Solve the system:
v = U \ ( L\ (s.*b) );

end
