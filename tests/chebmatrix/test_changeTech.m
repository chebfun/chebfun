function pass = test_changeTech(pref)
% Test CHEBMATRIX/CHANGETECH.

if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e-14;

%% CHEBTECH TO TRIGTECH.
f1 = chebfun(@(x) cos(x), [0 2*pi]);
f2 = chebfun(@(x) sin(x), [0 2*pi]); 
F = [ f1; f2 ];
G = changeTech(F, @trigtech);
pass(1) = isequal(get(G{1}.funs{1}, 'tech'), @trigtech);
pass(2) = isequal(get(G{2}.funs{1}, 'tech'), @trigtech);
pass(3) = norm(F{1} - G{1}, inf) < tol;
pass(4) = norm(F{2} - G{2}, inf) < tol;

%% TRIGTECH TO CHEBTECH.
f1 = chebfun(@(x) cos(x), [0 2*pi], 'trig');
f2 = chebfun(@(x) sin(x), [0 2*pi], 'trig'); 
F = [ f1; f2 ];
tech = pref.tech;
G = changeTech(F, tech);
pass(5) = isequal(get(G{1}.funs{1}, 'tech'), tech);
pass(6) = isequal(get(G{2}.funs{1}, 'tech'), tech);
pass(7) = norm(F{1} - G{1}, inf) < tol;
pass(8) = norm(F{2} - G{2}, inf) < tol;

%% TEST WITH MIXED OBJETCS.
f1 = 1;
f2 = chebfun(@(x) cos(10*x), [0 2*pi], 'trig');
f3 = chebfun(@(x) sin(10*x), [0 2*pi]);
F = [ f1; f2; f3 ];
G = changeTech(F, @trigtech);
% A scalar must stay a scalar:
pass(9) = isnumeric(G{1});
pass(10) = isequal(G{1}, 1);
% F{2} was already of the right tech, so we didn't do anything:
pass(11) = ( norm(F{2} - G{2}, inf ) == 0 );
% F{3} should have been converted:
pass(12) = isequal(get(G{3}.funs{1}, 'tech'), @trigtech);
pass(13) = norm(F{3} - G{3}, inf) < tol;

end