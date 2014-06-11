function bc = createBC(bcArg, ends)
%CREATEBC converts boundary condition syntax to chebfuns so constructBC.m
% can deal with them.
% 
% This command attempts to convert the "boundary" data into chebfuns. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isa(bcArg,'function_handle') )
    if ( nargin(bcArg) <= 1 )
        % convert to a chebfun.
        bc = chebfun( bcArg , ends );
    elseif ( nargin(bcArg) == 2 )
        % This allows for N.lbc = @(x,u) diff(u) - x + 1. 
        % We cannot convert this to a CHEBFUN so pass it on to MLDIVIDE.
        bc = bcArg;
    end
elseif ( isa(bcArg,'double') )
    % This allows for N.lbc = DOUBLE. 
    bc = chebfun( bcArg, ends );
elseif ( isa(bcArg,'chebfun') )
    % This allows for N.lbc = CHEBFUN. 
    bc = bcArg;
elseif ( isa(bcArg,'char') )
    if ( strcmpi(bcArg, 'periodic') )
        % Pass on to MLDIVIDE
        bc = 'periodic'; 
    else
        error('Chebop2:bcArg:word','Unrecognised bc string');
    end
else
    error('Chebop2:bcArg:type','Unrecognised boundary condition syntax');
end

end