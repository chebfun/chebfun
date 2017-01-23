function display(F)
%DISPLAY   Display a CHEBFUN3V object.
%   DISPLAY(F) outputs important information about the CHEBFUN3V object F 
%   to the command line, including its domain of definition and a summary 
%   of its structure.
%
%   This is called automatically when the semicolon is not used at the end 
%   of a statement that results in a CHEBFUN3V object.
%
% See also CHEBFUN3V/DISP and CHEBFUN3/DISPLAY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isequal(get(0, 'FormatSpacing'), 'compact') )
	disp([inputname(1), ' =']);
	disp(F);
else
	disp(' ');
	disp([inputname(1), ' =']);
	disp(' ');
    disp(F);
end

end