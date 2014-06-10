function M = feval(L, n, flag)
%FEVAL   Deprecated function, provided for limited backward compatability.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

warning('chebfun:linop:fevalDeprecated',...
    ['This function is provided only for limited backward compatibility.',...
    ' Use MATRIX instead.'] );
warning('off', 'chebfun:linop:fevalDeprecated')

if ( prod(size(L)) > 1 )
    error('This syntax is not available for multivariable systems.')
end

if ( nargin < 3 )
    flag = 'bc';
end

if ( isa(n, 'chebfun') )
    M = L*n;   % application to a function
else
    % Use a colloc2 discretization.
    disc = colloc2(L);
    disc.dimension = n;
    [PA, P, B, A] = matrix(disc);
    
    % Depending on the flag, we will do different things about boundary
    % conditions.
    switch(flag)
        case 'nobc'
            M = A;
        case 'bc'
            M = [B; PA];
        otherwise  % oldschool

            k = size(B,1);     % number of rows to drop
            k2 = ceil(k/2);    % about half for top
            try 
                A = cell2mat(A);
            catch
                error('oldschool does not support this problem');
            end

            A(1:k2,:) = B(1:k2,:);
            krem = k - k2;          % number remaining
            A(end-krem+1:end, :) = B(k2+1:end, :);
            M = A;
    end
end
        
end

