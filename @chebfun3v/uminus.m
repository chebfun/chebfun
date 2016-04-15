function F = uminus( F )
%- Unary minus of a CHEBFUN3V
%   -F returns the CHEBFUN3V negated componentwise. 
%
%   UMINUS(F) is called by the syntax -F. 

% Empty check: 
if (isempty( F ) )
    return
end

% Take uminus of each component:
F.components = cellfun( @uminus, F.components, 'UniformOutput', false );

end