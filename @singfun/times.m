function s = times(f,g)
%.* Multiply two singfuns

% This method will be called only if both F and G are SINGFUNS or at the most
% one of F and G is a scalar double.

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
    s = g;
    s.smoothPart = f * g.smoothPart;
    return
end

if ( isa(g,'double') )
    s = f;
    s.smoothPart = g * f.smoothPart;
    return
end

%%
% mutliplication of two singfuns
s = singfun;
% multiply the smooth parts
s.smoothPart = (f.smoothPart).*(g.smoothPart);
% add the exponents
s.exponents = f.exponents + g.exponents;

%%
% Check if after multiplication the 
% type of singularity has changed
% or if it can be removed.
% [TODO]: Since exponents are negative,
% it's impossible to remove a singularity
% by adding the exponents?
tol = singfun.pref.singfun.eps;
if ( abs(s.exponents(1)) > 100*tol )
    s.isSingEnd(1) = 1;
    if ( abs(s.exponents(1) - round(s.exponents)) < 100*tol )
        s.singType{1} = 'pole';
    else
        s.singType{1} = 'branch';
    end
else
    s.isSingEnd(1) = 0;
    s.singType{1} = 'none';
end

if ( abs(s.exponents(2)) > 100*tol )
    s.isSingEnd(2) = 1;
    if ( abs(s.exponents(2) - round(s.exponents)) < 100*tol )
        s.singType{2} = 'pole';
    else
        s.singType{2} = 'branch';
    end
else
    s.isSingEnd(2) = 0;
    s.singType{2} = 'none';
end
