function disp(F)
%DISP   Display a CHEBFUN3V.
%
% See also DISPLAY.

loose = strcmp(get(0,'FormatSpacing'),'loose');

% Compact version
if ( isempty( F ) )
    fprintf('empty chebfun3v\n')
    return
end

if ( F.isTransposed )
    tString = 'Row vector';
else
    tString = 'Column vector';
end

disp(['   chebfun3v object ' '(' tString ') containing:' ])
if ( loose )
    fprintf('\n');
end

% Display its CHEBFUN3 parts.
for j = 1:F.nComponents
    disp( F.components{j} );
end

end