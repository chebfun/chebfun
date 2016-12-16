function [err, lin] = valueTesting(f, numOut)
%VALUETESTING   Test that ADCHEBFUN is calling the correct method for the
%   function part of the methods.
%
%   ERR = VALUETESTING(F, NUMOUT) tests that ADCHEBFUN is calling the correct
%   method for the function part of the methods.
%
%   Here, the inputs are:
%   F       -- a function handle
%   NUMOUT  -- an optional argument, used for methods with more than one output
%              (in particular, ellipj)
%
%   and the output is
%   ERR     --  the infinity norm of the difference between applying F to a
%               CHEBFUN and an ADCHEBFUN.
%   LIN     --  linearity information. LIN == 1 if the method is determined to
%               be linear, 0 otherwise.
%
% See also: TAYLORTESTING, TAYLORTESTINGBINARY, VALUETESTINGBINARY.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This method proceeds as follows:
%   1. Construct an arbitrary CHEBFUN U, and a corresponding ADCHEBFUN V.
%   2. Evaluate the function handle F on both U and V.
%   3. Return the infinity norm of the difference between the results (which we
%      should expect to be zero).
%   4. Also return the variable LIN, that contains linearity information for the
%      function handle we were evaluating.

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

% Is the output linear?
lin = cellfun(@(fv) isLinear(fv), fv);

end
