%% Encryption of a message with SCRIBBLE
% Nick Trefethen, 10th April 2012

%%
% (Chebfun example fun/Encryption.m)
% [Tags: #encryption, #message]

%%
% SCRIBBLE produces piecewise linear complex chebfun whose plots look like
% words, like this:
message = scribble('This is the message');
LW = 'linewidth'; lw = 2;
plot(message, LW, lw), axis equal

%%
% Here's another string:
key = scribble('Aardvarks eat ants');
plot(key, 'r', LW, lw), axis equal

%%
% Now if we plot the sum of the two, we get nonsense:
encrypted = message + key;
plot(encrypted, 'm', LW, lw), axis equal

%%
% So we've invented a new encryption scheme!  For of course the original message
% can be recovered by subtracting off that key:
message2 = encrypted - key;
plot(message2, LW, lw), axis equal

%%
% So long as we're investigating the world's most expensive and least secure
% method of encryption, we might as well tangle up the text in the complex plane
% a bit too.  I'll bet you can't read this:
scrambled = exp(1.5i*(encrypted));
plot(scrambled, 'g', LW, lw), axis equal

%%
% But we can get the message back with a little unscramble:
message3 = unwrap(log(scrambled))/1.5i - 1 - key;
plot(message3, LW, lw), axis equal
