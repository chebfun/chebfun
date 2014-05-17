function splitting(on_off)
%SPLITTING   CHEBFUN splitting option
%   SPLITTING ON allows the Chebfun constructor to split the interval by a
%   process of automatic subdivision and edge detection.  This option is
%   recommended when working with functions with singularities. 

%   SPLITTING OFF disables this kind of automatic splitting, and is
%   recommended for working with functions that are complicated but still
%   smooth. Even with splitting off, breakpoints may still be introduced by
%   the MAX, MIN, ABS, CEIL, FLOOR, and ROUND commands.  One may switch
%   freely back and forth between the two modes during a Chebfun
%   computation. 

%   SPLITTING is OFF by default, by itself, displays the current splitting
%   state,

%   Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if nargin==0 
    switch chebfunpref('splitting')
        case 1 
            disp('SPLITTING is currently ON')
        case 0
            disp('SPLITTING is currently OFF')
    end
else
    if strcmpi(on_off, 'on')
        chebfunpref('splitting',true)
    elseif strcmpi(on_off, 'off') 
        chebfunpref('splitting',false)
    else
        error('CHEBFUN:split:UnknownOption',...
          'Unknown splitting option: only ON and OFF are valid options.')
    end
end
