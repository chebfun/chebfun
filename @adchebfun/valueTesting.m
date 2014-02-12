function err = valueTesting(f, numOut)
%VALUETESTING  Test that ADCHEBFUN is calling the correct method for the
%   function part of the methods.
%
% Here:
%   F is a function handle
%   numOut is an optional argument, used for methods with more than one outputs
%       (in particular, ellipj)

% TODO: Describe algorithm.

% Default value of NUMOUT
if ( nargin == 1)
    numOut = 1;
end

% Seed random generator to ensure same values.
seedRNG(6179);

% Length of test functions:
N = 8;

% Generate an arbitrary CHEBFUN to evaluate the function at:
u = 0.1*chebfun(rand(N, 1)) + .5;

% Construct a corresponding ADCHEBFUN:
v = adchebfun(u);

% Evaluate f at both u and v. FU will be a CHEBFUN, while FV will be an
% ADCHEBFUN
[fu{1:numOut}] = f(u);
[fv{1:numOut}] = f(v);

% We should expect that the FUNC field of FV matches the FU
err = max(cellfun(@(fu, fv) norm(fu - fv.func,'inf'), fu, fv));

end