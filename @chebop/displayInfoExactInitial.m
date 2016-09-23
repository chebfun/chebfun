function displayInfoExactInitial(pref)
%DISPLAYINFOEXACTINITIAL    Called if we pass CHEBOP/SOLVEBVP an initial guess
%                           that solves the BVP.
%
% [displayInfoExactInitial(U0, PREF) prints out information before the Newton
% iteration starts. This method is called if the initial guess passed to
% CHEBOP/SOLVEBVP appears to be a solution to the problem. In particular if
% PREF.display = 'iter', it prints the information that the initial guess
% solved the BVP
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Obtain preferences for what we want to show
display  = pref.display;

% Do we want to print to the command window?
if ( strcmpi(display, 'iter') )

    % Create a long string of dashes...
    dashString = repmat('-', 1, 62);
    % ... and print to the command window
    fprintf('\n%s\n', dashString);
    
    % Print to the command window:
    fprintf('Initial guess appears to be a solution of the BVP. \n');
    fprintf('Initial guess returned as the solution.\n');

    % Print the long string of dashes again:
    fprintf('%s\n', dashString);
end

end
