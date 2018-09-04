function disp(f)
%DISP  Display a BALLFUNV to the command line.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

disp('   ballfunv object containing')
fprintf('\n');
F = f.comp;
for j = 1:3
    disp(F{j});
end
end