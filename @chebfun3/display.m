function display(F)
%DISPLAY    Display information about a CHEBFUN3.
%   DISPLAY(F) outputs information about the CHEBFUN3 F to the 
%   command window, including its domain of definition, trilinear rank, 
%   and a summary of its structure.
%
%   It is called automatically when the semicolon is not used at the end of
%   a statement that results in a CHEBFUN3.
%
% See also CHEBFUN3/DISP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
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