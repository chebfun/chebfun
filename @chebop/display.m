function display(X)
%DISPLAY   Display information about a CHEBOP.
%   DISPLAY(F) outputs important information about the CHEBOP F to the command
%   window. DISPLAY(F) is called automatically when the semicolon is not used at
%   the end of a statement that results in a CHEBOP.
%
% See also DISP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
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
