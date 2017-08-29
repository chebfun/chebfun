function disp(F)
%DISP   Display a DISKFUN to the command line.
% 
% See also DISKFUN/DISPLAY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

loose = strcmp( get(0, 'FormatSpacing'), 'loose' );

% Get display style and remove trivial empty DISKFUN case. 
if ( isempty(F) )
    fprintf('    empty diskfun\n')
    if ( loose )
        fprintf('\n');
    end
    return
end

% Get information that we want to display:
len = length(F);                          % Numerical rank
vscl = vscale(F);                         % vertical scale

% Display the information: 
disp('     diskfun object ')
fprintf('       domain        rank    vertical scale\n');
fprintf('      unit disk   %6i          %3.2g\n', len, vscl);

if ( loose )
    fprintf('\n');
end

end
