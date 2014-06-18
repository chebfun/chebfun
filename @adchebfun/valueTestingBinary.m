function [err, lin] = valueTestingBinary(func)
%ERR = VALUETESTING(F)   Test that ADCHEBFUN is calling the correct method for
%   the function part of the methods. This method is used for binary operators.
%
% Here, the input is:
%   F       -- a function handle
%
% and the output is
%   ERR     --  a vector containing the infinity norm of the difference between
%               applying F to CHEBFUNS and an ADCHEBFUNS.
%   LIN     --  linearity information. LIN(i) == 1 if the method is determined
%               to be linear for the ith combination of variables, 0 otherwise.
%          
%
% See also: TAYLORTESTING, TAYLORTESTINGBINARY, VALUETESTING.

% This method proceeds as follows:
%   1. Construct arbitrary CHEBFUNS U1 and U2, and corresponding ADCHEBFUNS V1
%      and V2. Also construct two other arbitrary CHEBFUNS W1 and W2, as well as
%      arbitrary scalars S1 and S2.
%   2. Evaluate the function handle F on various combinations of U1 and U2 along
%      with W1, W2, S1 and S2, and then do the same operations for V1 and V2. 
%   3. Return a vector containing the infinity norm of the difference between
%      the results of the matching operations (which we should expect to be zero
%      in all cases).

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


%% Initialise

% Seed random generator to ensure same values.
seedRNG(6179);

% Length of test functions:
N = 8;

% Generate two arbitrary CHEBFUN objects to evaluate the function at:
u1 = chebfun(0.1*rand(N, 1) + .5);
u2 = chebfun(0.1*rand(N,1) + .5);

% Construct corresponding ADCHEBFUN objects. Here, we seed the ADCHEBFUN objects
% to ensure that the dimensions of their derivatives are correct:
v1 = seed(adchebfun(u1), 1, 2);
v2 = seed(adchebfun(u2), 2, 2);

% Also create two arbitrary CHEBFUNS to test with:
w1 = chebfun(0.1*rand(N, 1) + .5);
w2 = chebfun(0.1*rand(N,1) + .5);

% And finally two scalars to test with:
s1 = rand();
s2 = rand();

% Initialise error and linearity vectors:
err = zeros(1,5);
lin = err;
%% Create various combinations

% ADCHEBFUN and ADCHEBFUN:
v1v2 = func(v1, v2);

% ADCHEBFUN and CHEBFUN:
v1w2 = func(v1, w2);
w1v2 = func(w1, v2);

% ADCHEBFUN and SCALAR:
v1s2 = func(v1, s2);
s1v2 = func(s1, v2);

% Do CHEBFUN counterparts for allowing comparing values:
u1u2 = func(u1, u2);
u1w2 = func(u1, w2);
w1u2 = func(w1, u2);
u1s2 = func(u1, s2);
s1u2 = func(s1, u2);

%% Compare values
err(1) = norm(v1v2.func - u1u2);
err(2) = norm(v1w2.func - u1w2);
err(3) = norm(w1v2.func - w1u2);
err(4) = norm(v1s2.func - u1s2);
err(5) = norm(s1v2.func - s1u2);

%% Linearity information
lin(1) = isLinear(v1v2);
lin(2) = isLinear(v1w2);
lin(3) = isLinear(w1v2);
lin(4) = isLinear(v1s2);
lin(5) = isLinear(s1v2);
end
