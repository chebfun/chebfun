function pass = test_twocomponents(pref)
% Testing chebfun3v objects with two components. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;

% Different calls to the constructor. 
f = chebfun3(@(x,y,z) x); 
F1 = chebfun3v(f, f); % from chebfun3 objects. 
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x); % from handles.
pass(1) = norm(F1-F2) < tol;

% With domains. 
dom = [-1 1 -1 1 -1 1];
f = chebfun3(@(x,y,z) x, dom);
F1 = chebfun3v(f, f, dom); % from chebfun3 objects. 
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x, dom); % from handles.
pass(2) = norm(F1-F2) < tol;

% Plus chebfun3v objects. 
G = chebfun3v(f, 2*f);
H = F1 + G;
pass(3) = norm(H(1, 1, 1) - [2 3]')  < tol;

H = G + 1;
pass(4) = norm( H(pi/6, 1, 1) - pi/6*[1 2]' - 1 ) < tol;

H = G + [1 2];
pass(5) = norm(H(pi/6, 1, 1) - pi/6*[1 2]' -[1 2]') < tol;

% in reverse
H = 1 + G;
pass(6) = norm(H(pi/6, 1, 1) - pi/6*[1 2]' - 1) < tol;

H = [1 2] + G; 
pass(7) = norm(H(pi/6, 1, 1) - pi/6*[1 2]' -[1 2]') < tol;

% .* chebfun3v objects. 
G = chebfun3v(f, 2*f);
H = G.*1;
pass(8) = norm(H(pi/6, 1, 1) - G(pi/6, 1, 1)) < tol;

H = G.*[1 2]'; 
pass(9) = norm(H(pi/6, 1, 1) -  pi/6*[1 4]') < tol;

% in reverse
H = 1.*G; 
pass(10) = norm( H(pi/6, 1, 1) - G(pi/6, 1, 1) ) < tol;

H = [1 2]'.*G; 
pass(11) = norm(H(pi/6, 1, 1) -  pi/6*[1 4]') < tol;

% * chebfun3v objects. 
G = chebfun3v(f, 2*f);
H = G*1;
pass(12) = norm(H(pi/6, 1, 1) - G(pi/6, 1, 1)) < tol;

H = G'*[1 2]';   % inner product
pass(13) = norm(H - (f+4*f)) <tol;

% in reverse
H = 1*G;
pass(14) = norm(H(pi/6, 1, 1) - G(pi/6, 1, 1)) < tol;

try 
    H = [1 2]'*G;  % this should fail .
catch
    pass(15) = 1;
end

end