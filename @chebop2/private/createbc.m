function bc = createbc(bcArg, ends)
%CREATEBC converts boundary condition syntax to chebfuns so constructbc.m
% can deal with them.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isa(bcArg,'function_handle') )
    if ( nargin(bcArg) <= 1 )
        % convert to a chebfun.
        bc = chebfun( bcArg , ends );
    elseif ( nargin(bcArg) == 2 )
        % such as @(x,u) diff(u) - x + 1, leave until mldivide...
        bc = bcArg;
    end
elseif ( isa(bcArg,'double') )
    bc = chebfun( bcArg, ends );
elseif ( isa(bcArg,'chebfun') )
    bc = bcArg;
elseif ( isa(bcArg,'char') )
    if ( strcmpi(bcArg, 'periodic') )
        bc = 'periodic'; 
    else
        error('Chebop2:bcArg:word','Unrecognised bc string');
    end
else
    error('Chebop2:bcArg:type','Unrecognised boundary condition syntax');
end

end