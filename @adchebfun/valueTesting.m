function err = valueTesting(f, numOut)
% VALUETESTING  Test that ADCHEBFUN is calling the correct method for the
%   function part of the methods.
% Here:
%   F is a function handle
%   numOut is an optional argument, used for methods with more than one outputs
%       (in particular, ellipj)

% Default value of NUMOUT
if ( nargin == 1)
    numOut = 1;
end

% Seed random generator to ensure same values.
seedRNG(6179);

% Generate an arbitrary CHEBFUN to evaluate the function at
u = 0.1*chebfun(rand(15,1)) + .5;

% Construct a corresponding ADCHEBFUN
v = adchebfun(u);

% Call the method
if ( numOut == 1)
    
    % Evaluate f at both u and v. FU will be a CHEBFUN, while FV will be an
    % ADCHEBFUN
    fu = f(u);
    fv = f(v);
    
    % We should expect that the FUNC field of FV matches the FU
    err = norm(fu-fv.func,'inf');
    
elseif ( numOut == 3)
    % Evaluate f at both u and v. FU, GU, HU will be CHEBFUN objects, while FV,
    % GV and HV will be ADCHEBFUN objects
    [fu, gu, hu] = f(u);
    [fv, gv, hv] = f(v);
    
    % We should expect that the FUNC fields of *V matches those of *U
    err = max([norm(fu-fv.func, 'inf'), norm(gu-gv.func, 'inf'), ...
        norm(hu-hv.func, 'inf')]);
else
    error('CHEBFUN:ADCHEBFUN:valueTesting', ...
        'Unexpected number of output arguments')
end