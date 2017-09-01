function out = isempty( F )
%ISEMPTY empty boolean check for a DISKFUNV object. 
%   ISEMPTY(F) returns 1 if every component of F is an empty DISKFUN, and
%   return 0 otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Definite an empty DISKFUNV if it's components field is empty: 
if ( isempty( F.components ) ) 
    out = 1; 
    return
end

% Also an empty DISKFUNV if the components field is not empty, but contains
% empty fields. We check this by calling isempty of each component:
out = cellfun( @isempty, F.components, 'UniformOutput', false );
out = all(cell2mat( out ) );

end