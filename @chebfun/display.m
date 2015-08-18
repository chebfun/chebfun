function display(X)
%DISPLAY   Display information about a CHEBFUN.
%   DISPLAY(F) outputs important information about the CHEBFUN F to the command
%   window, including its domain of definition, its length (number of sample
%   values used to represent it), and a summary of its values at its endpoints.
%   DISPLAY(F) is called automatically when the semicolon is not used at the end
%   of a statement that results in a CHEBFUN.
%
%   Note that the VSCALE displayed is that returned by VSCALE(F, 'SUP').
%
% See also DISP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isequal(get(0, 'FormatSpacing'), 'compact') )
    disp([inputname(1), ' =']);
    disp(X);
else
    fprintf('\n');
    disp([inputname(1), ' =']);
    fprintf('\n');
    disp(X);
end

end
