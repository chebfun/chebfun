function s = rdivide(f, g)
%./   Divide SINGFUNS with SINGFUNS
%
%   This method will be called only if both F and G are SINGFUNS or at the 
%   most one of F and G is a scalar double.
%
% See also LDIVIDE, TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Case of empty arguments.
if ( isempty(f) || isempty(g) )
    % Return an empty SINGFUN:
    s = singfun;
    return;
end

% Check if inputs are other than SINGFUNS or doubles.
if ( (~isa(f, 'singfun') && ~isa(f, 'double')) || ...
     (~isa(g, 'singfun') && ~isa(g, 'double')) )   
    error('SINGFUN:rdivide:Input can only be a singfun or a double')
end

%% Trivial Case:
% Scalar division by a double.
if ( isa(g,'double') )
    % Copy the other input (a SINGUN) in the output:
    s = f;
    % Divide the smooth part with the double g and return:
    s.smoothPart = (1/g) * f.smoothPart;
    return
end
%% Reciprocal:
% Reciprocal of a SINGFUN scaled by the double F.
if ( isa(f,'double') )    
    % Convert f to a SINGFUN and call rdivide again.
    
    % Make a zero SINGFUN:
    temp = singfun.zeroSingFun();        
    % Assign f as it's smooth part
    temp.smoothPart = singfun.constructSmoothPart(f, f, 1, []);
    % Change f to a SINGFUN:
    f = temp;
    % Call SINGFUN.RDIVIDE again
    s = f./g;
    return
end



%% SINGFUN./SINGFUN
% Check if g has any roots in the open interval (-1, 1)
r = roots( g.smoothPart );
% remove roots at the end points.
r = setdiff(r, [-1,1]);
if ( ~isempty(r) )
    error('SINGFUN:rdivide:Divide by zero error')
end

% Note: Exponents of f and g can all be zero to generate a singular function.
% Example: f = 1; g = cos(pi/2*x) with trivial exponents. Then s = f./g is 
% singular with non trivial exponents. So the result of f./g in general is a
% generic SINGFUN with possibly non-trivial exponents.

% construct the SINGFUN by a direct call to the constructor:
s = singfun(@(x) feval(f, x)./feval(g, x));

end