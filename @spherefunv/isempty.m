function out = isempty( F )
%ISEMPTY empty boolean check for a SPHEREFUNV object. 
%   ISEMPTY(F) returns 1 if every component of F is an empty SPHEREFUN, and
%   return 0 otherwise.

if ( isempty( F.components ) ) 
    out = 1; 
    return
end

% Take isempty of each component:
out = cellfun( @isempty, F.components, 'UniformOutput', false );
out = all(cell2mat( out ) );

end
