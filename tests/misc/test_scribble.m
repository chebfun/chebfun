% Test file for scribble.m.

function pass = test_scribble(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Save the current warning state and then clear it in case scribble() issued a
% warning previously.
[lastmsg, lastid] = lastwarn();
lastwarn('');

% To pass, this call to scribble() must not crash or generate any warnings.
c = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.,;:?!''"`_ )([]{}-+*/^=<>\|%#~@$&';
s = scribble(c);

[wmsg, wid] = lastwarn();
if ( strcmp(wid, 'CHEBFUN:scribble:unknownchar') )
    pass(1) = false;
else
    % If this run of scribble() didn't issue a warning, pass and restore the
    % warning state.
    pass(1) = true;
    lastwarn(lastmsg, lastid);
end

end
