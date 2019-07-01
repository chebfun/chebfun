function disp( f )
%DISP Display a BALLFUN object to the command line.
% 
% See also DISPLAY.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Nfunctions = numel(f);
% Vectorize the array
f = f(:);

% Display the object
if ( Nfunctions == 1 )
    % Get information that we want to display:
    display_ballfun(f);
else
    disp('   ballfun array containing')
    fprintf('\n');
    for i = 1:Nfunctions
        display_ballfun(f(i));
    end
end

end

% Display one ballfun function 
function display_ballfun(f)

loose = strcmp(get(0, 'FormatSpacing'), 'loose');

if isempty(f)
    fprintf('    empty ballfun\n')    
else
S = size(f);
vscl = vscale(f);
% Display the information:
disp('   ballfun object:')
fprintf('      domain           r    lambda    theta    vertical scale\n');
fprintf('     unit ball    %6i  %6i   %6i           %3.2g\n', S(1), S(2), S(3),vscl);
end

if ( loose )
    fprintf('\n');
end

end
