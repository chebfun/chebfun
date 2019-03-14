function disp(f)
%DISP  Display a BALLFUNV to the command line.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

loose = strcmp(get(0,'FormatSpacing'),'loose');

% Compact version:
if ( isempty(f) )
    fprintf('empty ballfunv\n\n')
    return
end

disp('   ballfunv object containing')
if ( loose )
    fprintf('\n');
end

F = f.comp;
for j = 1:3
    disp(F{j});
end

end