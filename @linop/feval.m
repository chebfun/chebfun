function M = feval(L, n, flag)
%FEVAL    Deprecated function, provided for limited backward compatibility.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:LINOP:feval:deprecated',...
    ['This function is provided only for limited backward compatibility.',...
    ' Use MATRIX instead.'] );
warning('off', 'CHEBFUN:LINOP:feval:deprecated')

if ( prod(size(L)) > 1 ) %#ok<PSIZE>
    error('CHEBFUN:LINOP:feval:multivariable', ...
        'This syntax is not available for multivariable systems.')
end

if ( nargin < 3 )
    flag = 'bc';
end
% TODO: Should we always use a chebcolloc2 representation here?
% if ( nargin < 4 )
%     pref = cheboppref();
%     discType = pref.discretization;
% elseif ( isa(discType, 'cheboppref') )
%     discType = discType.discretization;
% elseif ( ischar(discType) )
%     discType = str2func(discType);
% end
discType = @chebcolloc2;

if ( isa(n, 'chebfun') )
    M = L*n;   % application to a function
    
else

    disc = discType(L);
    disc.dimension = n;
    
    if ( strcmp(flag, 'oldschool') )
        disc.dimAdjust = 0;
    end
    
    [PA, ignored, B, A] = matrix(disc);
    
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
                error('CHEBFUN:LINOP:feval:oldschool', ...
                    'oldschool does not support this problem');
            end
            if ( k2 > 0 )
                A(1:k2,:) = B(1:k2,:);
                krem = k - k2;     % number remaining
                A(end-krem+1:end, :) = B(k2+1:end, :);
            end
            M = A;
    end
end
        
end

