function pass = test_twocomponents( pref ) 
% Testing chebfun3v objects with two components. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;
j = 1;

% Different calls to the constructor. 
f = chebfun3(@(x,y,z) x); 
F1 = chebfun3v(f, f); % from chebfun3 objects. 
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x); % from handles.
pass(j) = norm(F1-F2) < tol;
j = j+1; 

% With domains. 
dom = [-1 1 -1 1 -1 1];
f = chebfun3(@(x,y,z) x, dom);
F1 = chebfun3v(f, f, dom); % from chebfun3 objects. 
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x, dom); % from handles.
pass(j) = norm(F1-F2) < tol;
j = j+1;

% Plus chebfun3v objects. 
G = chebfun3v(f, 2*f);
H = F1 + G; 
pass(j) = norm(H(1, 1, 1) - [2 3]')  < tol;
j = j + 1; 

H = G + 1;
pass(j) = norm( H(pi/6, 1, 1) - pi/6*[1 2]' - 1 ) <tol;
j = j + 1;

H = G + [1 2];
pass(j) = norm(H(pi/6, 1, 1) - pi/6*[1 2]' -[1 2]') <tol;
j = j + 1;

% in reverse
H = 1 + G;
pass(j) = norm(H(pi/6, 1, 1) - pi/6*[1 2]' - 1) <tol; 
j = j + 1; 

H = [1 2] + G; 
pass(j) = norm(H(pi/6, 1, 1) - pi/6*[1 2]' -[1 2]') <tol; 
j = j + 1;

% .* chebfun3v objects. 
G = chebfun3v(f, 2*f);
H = G.*1;
pass(j) = norm(H(pi/6, 1, 1) - G(pi/6, 1, 1)) <tol; 
j = j + 1; 

H = G.*[1 2]'; 
pass(j) = norm(H(pi/6, 1, 1) -  pi/6*[1 4]') <tol;
j = j + 1;

% in reverse
H = 1.*G; 
pass(j) = norm( H(pi/6, 1, 1) - G(pi/6, 1, 1) ) <tol;
j = j + 1;

H = [1 2]'.*G; 
pass(j) = norm(H(pi/6, 1, 1) -  pi/6*[1 4]') <tol;
j = j + 1;

% * chebfun3v objects. 
G = chebfun3v(f, 2*f);
H = G*1; 
pass(j) = norm(H(pi/6, 1, 1) - G(pi/6, 1, 1)) <tol;
j = j + 1;

H = G'*[1 2]';   % inner product
pass(j) = norm(H - (f+4*f)) <tol;
j = j + 1;

% in reverse
H = 1*G;
pass(j) = norm(H(pi/6, 1, 1) - G(pi/6, 1, 1)) < tol;
j = j + 1;
try 
    H = [1 2]'*G;  % this should fail .
catch
    pass(j) = 1 ;
    j=j+1;
end

end