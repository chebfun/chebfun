function disp(f)
%DISP Display a BALLFUN function to the command line

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Nfunctions = numel(f);
% Vectorize the array
f = f(:);

% Display the object
if Nfunctions == 1
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
S = size(f);
% Display the information:
disp('   ballfun object:')
fprintf('      domain           r    lambda    theta\n');
fprintf('     unit ball    %6i    %6i   %6i\n', S(1), S(2), S(3));
end
