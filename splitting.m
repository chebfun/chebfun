function varargout = splitting(on_off)
%SPLITTING   CHEBFUN splitting option.
%   SPLITTING('ON') allows the CHEBFUN constructor to split the interval by a
%   process of automatic subdivision and edge detection. This option is
%   recommended when working with functions with singularities or jumps.
%
%   SPLITTING('OFF') disables this kind of automatic splitting, and is
%   recommended for working with functions that are complicated but still
%   smooth. Even with splitting off, breakpoints may still be introduced by the
%   MAX, MIN, ABS, CEIL, FLOOR, and ROUND commands. One may switch freely back
%   and forth between the two modes during a Chebfun computation.
%
%   SPLITSTATE = SPLITTING(...) will return the current state splitting option
%   (as a string), before any changes are applied.
%
%   SPLITTING by itself, displays a string detailing the current splitting
%   state. By default, this is OFF.
%
% See also CHEBFUNPREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargout > 0 && nargin == 0 )
    % Return current splitting state:
    splitState = chebfunpref().splitting;
    
elseif ( nargin == 0 )
    % Display splitting state:
    switch ( chebfunpref().splitting )
        case 1
            disp('SPLITTING is currently ON.')
        case 0
            disp('SPLITTING is currently OFF.')
    end
    
else
    % Throw a warning:
    warning('CHEBFUN:splitting:deprecated', ...
        ['The syntax ''splitting on'' is deprecated.\n', ...
        'Please see CHEBFUNPREF documentation for further details.']);
    % But only throw it once:
    warning('off', 'CHEBFUN:splitting:deprecated');
    
    if ( nargout > 0 )
        splitState = chebfunpref().splitting;
    end
    
    if ( strcmpi(on_off, 'on') )
        chebfunpref.setDefaults('splitting', true);
        
    elseif strcmpi(on_off, 'off')
        chebfunpref.setDefaults('splitting', false);
        
    else
        error('CHEBFUN:splitting:UnknownOption',...
            'Unknown splitting option: only ON and OFF are valid options.')
    end
    
end

if ( nargout > 0 )
    onOffStr = {'off', 'on'};
    varargout{1} = onOffStr{splitState+1};
end

end
