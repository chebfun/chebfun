function bc = createBC(bcArg, ends)
%CREATEBC converts boundary condition syntax to chebfuns so constructBC.m
% can deal with them.
% 
%   BC = CREATEBC(BCARG, ENDS) usually constructs a chebfun for the
%   homogeneous part of the boundary conditions.  If the linear constraint
%   is sufficiently compliciated then this command gives up and passes
%   resposibility to the solver. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%%% Developer note %%%%
% This command attempts to convert the "boundary" data into chebfuns. It is
% called when the boundary conditions are first assigned to the chebop2 
% object. 

if ( isa(bcArg,'function_handle') )
    if ( nargin(bcArg) <= 1 )
        % convert to a chebfun.
        bc = chebfun( bcArg , ends );
    elseif ( nargin(bcArg) == 2 )
        % This allows for N.lbc = @(x,u) diff(u) - x + 1. 
        % We cannot convert this to a CHEBFUN so pass it on to the SOLVER.
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
        % Pass on to SOLVER
        bc = 'periodic'; 
    else
        error('Chebop2:bcArg:word','Unrecognised bc string');
    end
else
    error('Chebop2:bcArg:type','Unrecognised boundary condition syntax');
end

end