function disp(f)
%DISP  Display a BALLFUNV to the command line.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Nfunctions = numel(f);
% Vectorize the array
f = f(:);

% Display the object
if Nfunctions == 1
    % Get information that we want to display:
    display_one_ballfunv(f);
else
    disp('   ballfunv array containing')
    fprintf('\n');
    for i = 1:Nfunctions
        display_ballfunv(f(i));
    end
end
end

% Display one ballfunvCart 
function display_one_ballfunv(f)
disp('   ballfunv object containing')
fprintf('\n');
F = f.comp;
for j = 1:3
    disp(F{j});
end
end

% Display a ballfunv f
function display_ballfunv(f)
S = size(f);
% Display the information:
disp('   ballfunv object:')
fprintf('      domain           r    lambda    theta\n');
fprintf('     unit ball    %6i    %6i   %6i\n', S(1), S(2), S(3));
end
