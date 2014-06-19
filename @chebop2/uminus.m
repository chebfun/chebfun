function N = uminus(N)
%UMINUS   Unitary minus for CHEBOP2 objects.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Negate all the variables coefficents: 
A = N.coeffs; 
if ( iscell(A) )
    for jj = 1:size(A, 1)
        for kk = 1:size(A, 2) 
            A{jj,kk} = -A{jj,kk};
        end
    end
else
    A = -A; 
end

% Do not negate the BCs:
if ( ~isempty(N.lbc) || ~isempty(N.rbc) || ~isempty(N.ubc) || ~isempty(N.dbc) )
    warning('CHEBFUN:CHEBOP2:unimus:BCs', ...
        'Operator has BCs. These were not negated.');
end

% Update the properties of the CHEBOP2:
N.coeffs = A; 
op = N.op; 
N.op = @(u) -op(u); 

end
