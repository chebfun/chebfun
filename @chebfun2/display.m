function display(F)
%DISPLAY   Display a CHEBFUN2 to the command line.
% 
% DISPLAY(F) outputs important information about the CHEBFUN2 F to the
% command window, including its domain of definition, length (number of 
% pivots used to represent it), and a summary of its structure. 
%
% It is called automatically when the semicolon is not used at the
% end of a statement that results in a CHEBFUN2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%
% Get display style and remove trivial empty CHEBFUN2 case. 

loose = strcmp( get(0, 'FormatSpacing'), 'loose' );
if ( loose )
    fprintf('\n%s = \n\n', inputname(1))
else
    fprintf('%s = \n', inputname(1))
end

% Compact version:
if ( isempty(F) )
    fprintf('    empty chebfun2\n\n')
    return
end

% Get information that we want to display:
dom = F.domain;                           % Domain
len = length(F);                          % Numerical rank
[xx, yy]=meshgrid(dom(1:2), dom(3:4));   
vals = feval(F, xx, yy ).';               % Corner values
vals = vals(:);
vscl = vscale(F);                         % vertical scale

% Display the information: 

disp('chebfun2 object: (1 smooth surface)')
fprintf('       domain                 rank       corner values\n');

if ( isreal(vals) )
    fprintf('[%4.2g,%4.2g] x [%4.2g,%4.2g]   %6i     [%4.2g %4.2g %4.2g %4.2g]\n', dom(1), dom(2), dom(3), dom(4), len, vals);
else
    fprintf('[%4.2g,%4.2g] x [%4.2g,%4.2g]   %6i     [  complex values  ]\n', dom(1), dom(2), dom(3) , dom(4), len);
end
fprintf('vertical scale = %3.2g \n', vscl)

end
