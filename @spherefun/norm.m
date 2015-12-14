function [ normF, normloc ] = norm( f, p )
%NORM   Norm of a SPHEREFUN
% 
%    SPHEREFUN(F) = sqrt(integral of abs(F)^2).
%    SPHEREFUN(F, 2) = largest singular value of F.
%    SPHEREFUN(F,'fro') is the same as NORM(F).
%    SPHEREFUN(F, 1) = NOT SUPPORTED.
%    SPHEREFUN(F, inf) = global maximum in absolute value.
%    SPHEREFUN(F, max) = global maximum in absolute value.
%    SPHEREFUN(F, min) = NOT SUPPORTED
%
% Furthermore, the inf norm for SPHEREFUN objects also returns a second output,
% giving a position where the max occurs.

if ( nargin == 1 ) 
    % Default to 2-norm.
    p = 2;
end

if ( isempty( f ) )  
    % Empty chebfun has norm 0.
    normF = [];
    
else
    switch ( p )  % Different cases on different norms.
        case 1
            error('CHEBFUN:SPHEREFUN:norm:norm', ...
                'SPHEREFUN does not support L1-norm, yet');
            
        case {2, 'fro'}  % Definite integral of f.^2
           
            s = svd( f ); 
            normF = sqrt( sum2( s.^2 ) );  
            
        case {inf, 'inf', 'max'}
             [Y, X] = minandmax2(f);
             [normF, idx] = max( abs( Y ) );
             normloc = X( idx, : );
            
        case {-inf, '-inf', 'min'}
            error('CHEBFUN:SPHEREFUN:norm:norm', ...
                'SPHEREFUN does not support this norm.');
            
%         case {'op', 'operator'}
%             [C, D, R] = cdr( f ); 
%             L = C * D * R; 
%             s = svd( L ); 
%             normF = s(1);
            
        otherwise
           % TODO:
            error 
             if ( isnumeric(p) && isreal(p) )
                if ( abs(round(p) - p) < eps )
                     p = round(p); f = f.^p;
                     if ( ~mod(p,2) )
                        normF = ( sum2( f ) ).^( 1/p );
                     else
                         error('CHEBFUN:SPHEREFUN:norm:norm', ...
                             'p-norm must have p even for now.');
                    end
                 else
                     error('CHEBFUN:SPHEREFUN:norm:norm', ...
                         'SPHEREFUN does not support this norm.');
                 end
             else
                 error('CHEBFUN:SPHEREFUN:norm:unknown', 'Unknown norm.');
             end
            
    end
end

end