function display(X)
%DISPLAY   Display a DISKFUNV.
%   DISPLAY(F) outputs important information about the DISKFUNV F to the
%   command window, including its domain of definition, length (number of pivots
%   used to represent it), and a summary of its structure.
%
%   It is called automatically when the semicolon is not used at the end of a
%   statement that results in a DISKFUNV.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isequal(get(0, 'FormatSpacing'), 'compact') )
	disp([inputname(1), ' =']);
	disp(X);
else
	disp(' ');
	disp([inputname(1), ' =']);
	disp(' ');
    disp(X);
end

end
