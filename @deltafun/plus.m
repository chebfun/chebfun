function s = plus(f, g)
%+   Addition of DELTAFUN objects.
%   F + G adds F and G, where F and G may be DELTAFUN objects or scalars.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

s = deltafun;
% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    % return with the empty DELTAFUN:
    return;
end

% One of the arguments i.e. f or g is necessarily a DELTAFUN object. Otherwise, 
% this overloaded plus would not have been called.

% If g is a double, f is a deltafun:
if ( isa(g, 'double') )
    s = f;
    s.funPart = f.funPart + g;
end

% Recursive call, by making the deltafun g as the first argument:
if ( isa(f, 'double') )
    s = g + f;
    return
end

if ( isa(f, 'deltafun') && isa(g, 'deltafun') )    
    s.funPart = f.funPart + g.funPart;
    % This is the main part:
    s.delta = f.delta;
end

%% %%%%%%%%%%%%% DELTAFUN + CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should be removed, since CHEBFUNs should be casted to DELTAFUNs first and
% then added
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive call: If g is a chebfun, upgrade it to a deltafun and call plus
% again.
if ( isa(f, 'deltafun') && isa(g, 'chebfun') )
    t = deltafun.zeroDeltaFun();
    t.funPart = g;
    s = f + t;
    return
end

%%
if ( issmooth(s) )
    % return a smooth object if the operation has removed all 
    % delta functions.
    s = s.funPart;
end


end