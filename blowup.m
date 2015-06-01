function varargout = blowup(on_off)
%BLOWUP   CHEBFUN blowup option.
%   CHEBFUN offers limited support for function which diverge to infinity on
%   their domain. The 'blowup' flag will determine whether such behaviour is
%   detected automatically. In particular;
%       BLOWUP(0) (or BLOWUP('off')): bounded functions only (default)
%       BLOWUP(1) (or BLOWUP('on')):  poles are permitted (integer order)
%       BLOWUP(2): blowup of arbitrary orders permitted (experimental)
%
%   STATE = BLOWUP(...) will return the current state of the blowup option (as a
%   double), before any changes are applied.
%
%   BLOWUP by itself, displays a string detailing the current blowup state. By
%   default, this is 0.
%
%   If the singular behaviour of the function is known in advance, it is usually
%   beneficial to specify this manually. This is achieved by defining 'exps'
%   (for "exponents") in the constructor. For example:
%       f = chebfun('1./x',[0 4], 'exps', [-1,0]);
%       f = chebfun('sin(10*pi*x) + 1./x', [-2 0 4], 'exps', [0 -1 0]);
%
%   Functions with noninteger blow up can also be represented to some extent,
%   but be warned that this functionality is still experimental and fragile.
%   For example:
%       f = chebfun('exp(x)./sqrt(1-x.^2)', 'exps', [-.5 -.5])
%   Nor is it required that the singular behavior result in a pole:
%       f = chebfun('exp(x).*(1+x).^pi.*(1-x).^sqrt(2)', 'exps', [pi sqrt(2)])
%   When the order of the blowup is unknown, one can pass a NaN and CHEBFUN will
%   attempt to compute the order automatically. (Currently it will only look for
%   integer poles, i.e. BLOWUP == 1).
%       f = chebfun('exp(x).*sqrt(1+x)./(1-x).^4', 'exps', [.5 NaN])
%
%   The SPLITTING ON functionality can also be used to detect such blowups.
%   Below are a number of ways to construct a chebfun representation of the
%   gamma function on the interval [-4 4] (in order of reliabilty)
%       f = chebfun('gamma(x)', [-4:0 4], 'exps', [-ones(1,9) 0])
%       f = chebfun('gamma(x)', [-4:0 4], 'exps', [-ones(1,5) 0])
%       f = chebfun('gamma(x)', [-4:0 4], 'blowup', 'on')
%       f = chebfun('gamma(x)', [-4 4], 'splitting', 'on', 'blowup', 'on')
%
% See also CHEBFUNPREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This is default behavior for "blowup on":
default = 1;

% Parse the inputs if they are a string:
if ( nargin > 0 && ischar(on_off) )
    check = strcmpi(on_off, {'on', 'off', '0', '1', '2'});
    if ( check(1) )
        on_off = default;
    elseif ( check(2) )
        on_off = 0;
    elseif ( any(check(3:end)) )        
        on_off = eval(on_off);
    else
        error('CHEBFUN:blowup:charIn', 'Unknown string input to BLOWUP().');
    end
end

if ( nargout > 0 && nargin == 0 )
    % Simply return the state:
    blowupState = getBlowupState();
    
elseif ( nargin == 0 )
    % Display the information:
    switch ( getBlowupState() )
        case 0
            disp('BLOWUP = 0: Bounded functions only')
        case 1 
            disp('BLOWUP = 1: Poles are permitted (integer order)')
        case 2 
            disp(['BLOWUP = 2: ', ...
            'Blowup of integer or non-integer orders permitted (experimental)'])
    end
    
else
    % Adjust the state:
    
    % Throw a warning:
    warning('CHEBFUN:blowup:deprecated', ...
        ['The BLOWUP() method is deprecated.\n', ...
        'Please see CHEBFUNPREF documentation for further details.']);
    % But only throw it once:
    warning('off', 'CHEBFUN:blowup:deprecated');
    
    % Return previous state:
    if ( nargout > 0 )
        blowupState = getBlowupState();
    end
    
    if ( on_off == 0 )
        chebfunpref.setDefaults('blowup', false);
        
    elseif ( on_off == 1 )
        chebfunpref.setDefaults('blowup', true);
        chebfunpref.setDefaults({'blowupPrefs', 'defaultSingType'}, 'pole')
        
    elseif ( on_off == 2 )
        chebfunpref.setDefaults('blowup', true);
        chebfunpref.setDefaults({'blowupPrefs', 'defaultSingType'}, 'sing')
        
    else
        error('CHEBFUN:blowup:unknownOption',...
          'Unknown blowup option: only ON, OFF, 0, 1, 2 are valid options.')
    end
    
end

% Return the previous state:
if ( nargout > 0 )
    varargout{1} = blowupState;
end

end

function state = getBlowupState()
    state = chebfunpref().blowup;
    if ( state )
        type = chebfunpref().blowupPrefs.defaultSingType;
        if ( any(strcmpi(type, 'sing')) )
            state = 2;
        end
    end
end

