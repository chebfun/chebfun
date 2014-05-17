function varargout = blowup(on_off)
%BLOWUP     CHEBFUN blowup option
% Chebfun offers limited support for function which diverge to infinity on
% their domain. 
%
% The 'blowup' flag will determine whether such behaviour is detected
% automatically. In particular;
%   BLOWUP=0 (or 'off'): bounded functions only (default)
%   BLOWUP=1 (or 'on' ): poles are permitted (integer order)
%   BLOWUP=2           : blowup of arbitrary orders permitted (experimental)
%
% DEFUALT = BLOWUP(NaN) will return 1 or 2, depending on whether blowups
% can only have integer or arbitrary order. This can be manually adjusted
% in this mfile.
%
% If the singular behaviour of the function is known in advance, it is
% usually beneficial to specify this manually. This is achieved by defining
% 'exps' (for "exponents") in the constructor. 
% For example:
%       f = chebfun('1./x',[0 4],'exps',[-1,0]);
%       f = chebfun('sin(10*pi*x) + 1./x',[-2 0 4],'exps',[0 -1 0]);
%
% Functions with noninteger blow up can also be represented to some extent,
% (but warned that this functionality is still experimental and fragile).
% For example:
%       f = chebfun('exp(x)./sqrt(1-x.^2)','exps',[-.5 -.5])
% Nor is it required that the singular behavior result in a pole:
%       f = chebfun('exp(x).*(1+x).^pi.*(1-x).^sqrt(2)','exps',[pi sqrt(2)])
% When the order of the blowup is unknown, one can pass a NaN and chebfun
% will attempt to compute the order automatically. (Currently it will only
% look for integer poles, i.e. BLOWUP == 1).
%       f = chebfun('exp(x).*sqrt(1+x)./(1-x).^4','exps',[.5 NaN])
%
% The SPLITTING ON functionality can also be used to detect such blowups.
% Below are a number of ways to construct a chebfun representation of the 
% gamma function on the interval [-4 4] (in order of reliabilty)
%       f = chebfun('gamma(x)',[-4:0 4],'exps',[-ones(1,9) 0])
%       f = chebfun('gamma(x)',[-4:0 4],'exps',[-ones(1,5) 0])
%       f = chebfun('gamma(x)',[-4:0 4],'blowup','on')
%       f = chebfun('gamma(x)',[-4 4],'splitting','on','blowup','on');

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% This is default behavior for "blowup on"
default = 1;

if ischar(on_off)
    check = strcmpi(on_off,{'on','off'});
    if check(1)
        on_off = default;
    elseif check(2)
        on_off = 0;
    else
        error('CHEBFUN:blowup:charin','Unkown string input to blowup.m');
    end
end

if nargin == 1 && isnan(on_off)
    if nargout == 1, varargout{1} = default;  end
    return
end

if nargin == 0 
    switch chebfunpref('blowup');
        case 0
            disp('BLOWUP=0: bounded functions only')
        case 1 
            disp('BLOWUP=1: poles are permitted (integer order)')
        case 2 
            disp('BLOWUP=2: blowup of integer or non-integer orders permitted (experimental)')
    end
else
    if strcmpi(on_off, 'on')
        chebfunpref('blowup',default)
    elseif strcmpi(on_off, 'off') || all(on_off==0);
        chebfunpref('blowup',0)
    elseif strcmpi(on_off, '1') || all(on_off==1);
        chebfunpref('blowup',1)
    elseif strcmpi(on_off, '2') || all(on_off==2);
        chebfunpref('blowup',2)      
    else
        error('CHEBFUN:blowup:UnknownOption',...
          'Unknown blowup option: only ON, OFF, 0, 1, 2, and NaN are valid options.')
    end
end

if nargout == 1
    varargout{1} = chebfunpref('blowup');
end
