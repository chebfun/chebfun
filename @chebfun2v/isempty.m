function bol = isempty( F )
%ISEMPTY empty boolean check for a chebfun2v object. 
% 
% ISEMPTY(F) returns 1 if every component of F is an empty chebfun2, and 
% return 0 otherwise. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 


if ( isempty( F.components ) ) 
    bol = 1; 
    return
end

% Take isempty of each component:
bol = cellfun( @isempty, F.components, 'UniformOutput', false );
bol = all (cell2mat( bol ) );

end