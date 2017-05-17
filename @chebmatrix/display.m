function display(X)
%DISPLAY   Display information about a CHEBMATRIX.
%   DISPLAY(L) prints the size of the CHEBMATRIX L and a list of the block class
%   types. If java is enabled, the class types for each block are hyperlinked,
%   and clicking them will call the display method for the item in that block.
%
% See also DISP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

name = inputname(1);

if ( isequal(get(0, 'FormatSpacing'), 'compact') )
	disp([name, ' =']);
	disp(X, name);
else
	disp(' ');
	disp([name, ' =']);
	disp(' ');
    disp(X, name);
end

end


