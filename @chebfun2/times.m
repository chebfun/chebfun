function h = times(f, g)
% .*   CHEBFUN2 times. 
%  

if ( isa(f, 'chebfun2') )    % CHEBFUN2 .* ??? 
    if ( isa(g, 'double') )  % CHEBFUN2 .* DOUBLE
        h = mtimes(f, g); 
    elseif ( isa( g, 'chebfun2') )
        %domain = domain_check(f, g);
        domain = [-1 1 -1 1];
        h = chebfun2(@(x, y) feval(f, x, y).*feval(g, x, y), domain);
    else
       error('CHEBFUN:MTIMES:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type %s ' ...
           'and %s.'], class(f), class(g));
    end
else
    h = times(g, f);     
end


end 