function [ normF ] = norm( f, p )
%NORM  
% For DISKFUN objects:
%    DISKFUN(F) = sqrt(integral of abs(F)^2).
%    DISKFUN(F, 2) = largest singular value of F.
%    DISKFUN(F,'fro') is the same as NORM(F).
%    DISKFUN(F, 1) = NOT IMPLEMENTED.
%    DISKFUN(F, inf) = global maximum in absolute value.
%    DISKFUN(F, max) = global maximum in absolute value.
%    DISKFUN(F, min) = NOT IMPLEMENTED


% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 ) 
    % Default to 2-norm.
    p = 2;
end

if ( isempty( f ) )  
    % Empty diskfun has norm 0.
    normF = [];
    
else
    switch ( p )  % Different cases on different norms.
        case 1
            error('CHEBFUN:DISKFUN:norm:norm', ...
                'DISKFUN does not support L1-norm, yet');
            
        case {2, 'fro'}  % Definite integral of f.^2
            % 
            s = svd( f ); 
            normF = sqrt( sum( s.^2 ) );  
            
         case {inf, 'inf', 'max'}
             [Y, X] = minandmax2(f);
             [normF, idx] = max( abs( Y ) );
             normloc = X( idx, : );
            
        case {-inf, '-inf', 'min'}
            error('CHEBFUN:DISKFUN:norm:norm', ...
                'DISKFUN does not support this norm.');
            
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
                        error('CHEBFUN:DISKFUN:norm:norm', ...
                            'p-norm must have p even for now.');
                    end
                else
                    error('CHEBFUN:DISKFUN:norm:norm', ...
                        'DISKFUN does not support this norm.');
                end
            else
                error('CHEBFUN:DISKFUN:norm:unknown', 'Unknown norm.');
            end
            
    end
end

end