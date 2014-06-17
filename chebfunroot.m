function out = chebfunroot()
%CHEBFUNROOT   Root directory of Chebfun installation.
%   S = CHEBFUNROOT() returns a string that is the name of the directory where
%   Chebfun is installed.

out = fileparts(which('chebfunroot'));

end