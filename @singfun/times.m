function s = times(f,g)
%.*   Multiply SINGFUNS with SINGFUNS
%
%   This method will be called only if both F and G are SINGFUNS or at the 
%   most one of F and G is a scalar double.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Check if inputs are other than SINGFUNS or doubles
if ( ~isa(f, 'singfun') && ~isa(f, 'double') )
    error( 'SINGFUN:times:Input can only be a singfun or a double' )
end

if ( ~isa(g, 'singfun') && ~isa(g, 'double') )
    error( 'SINGFUN:times:Input can only be singfun or a double' )
end

%%
% scalar multiplication cases
if ( isa(f,'double') )
    % copy the other input (a SINGUN) in the output
    s = g;
    % multiply the smooth part with the double and return
    s.smoothPart = f * g.smoothPart;
    return
end

if ( isa(g,'double') )
    % copy the other input (a SINGUN) in the output
    s = f;
    % multiply the smooth part with the double and return
    s.smoothPart = g * f.smoothPart;
    return
end

%%
% mutliplication of two SINGFUNS
s = singfun;
% multiply the smooth parts
s.smoothPart = (f.smoothPart).*(g.smoothPart);
% add the exponents
s.exponents = f.exponents + g.exponents;

%%
% Check if after multiplication the type of singularity has changed or if 
% it can be removed.
% [TODO]: Since exponents are negative, it's impossible to remove a
% singularity after mutiplying two SINGFUNS?
tol = singfun.pref.singfun.eps;
% loop through each end
for k = 1:2
    if ( s.exponents(k) < -100*tol )
        if ( abs(s.exponents(k) - round(s.exponents(k))) < 100*tol )
            s.singType{k} = 'pole';
        else
            s.singType{k} = 'sing';
        end
    else
        s.singType{k} = 'none';
    end
end

end