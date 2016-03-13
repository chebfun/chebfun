function disp(F)
%DISP   Display a DISKFUNV.
%
% See also DISPLAY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

loose = strcmp(get(0,'FormatSpacing'),'loose');

% Compact version
if ( isempty( F ) )
    fprintf('empty diskfunv\n\n')
    return
end

if ( F.isTransposed )
    tString = 'Row vector';
else
    tString = 'Column vector';
end

disp(['diskfunv object containing' ])
if ( loose )
    fprintf('\n');
end

% Display its two DISKFUN components.
for j = 1:F.nComponents
    disp( F.components{j} );
end

end
