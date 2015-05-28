function h = plus( f, g ) 
% PLUS      plus for SPHEREFUN objects. 
% 
% F + G adds F and G. F and G can be scalars or SPHEREFUN objects.

% Not implemented YET.  Need to do compress_plus.

if ( ~isa(f, 'spherefun') ) % ??? + SPHEREFUN
    
    h = plus(g, f);
    
elseif ( isempty(g) || isempty(f) ) % SPHEREFUN + []
    
    % Return empty SPHEREFUN.
    h = spherefun( );
    
elseif ( isa( g, 'double' ) )           % SPHEREFUN + DOUBLE
    
    % Convert g to a SPHEREFUN.
    g = spherefun( g );
    h = plus( f, g );
    
elseif ( ~isa(g, 'spherefun') )          % SPHEREFUN + ???
    
    error( 'SPHEREFUN:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class( f ), class( g ));
    
else                                     % SPHEREFUN + SPHEREFUN
      % Amalgamate the columns and rows of f and g
       h = f; 
       h.cols = [f.cols g.cols]; 
       h.rows = [f.rows g.rows];
       Bf = f.blockDiag; 
       Bg = g.blockDiag; 
       h.blockDiag = [Bf                    zeros(size(Bf,1),size(Bg,2)); 
                zeros(size(Bg,1), size(Bf,2))               Bg       ];
%        % Now compress: 
%        h = compress( h ); 
end

end