function varargout = scribble2( s, rk )
%SCRIBBLE2  Writing text with chebfun2 objects
%  SCRIBBLE('STRING') returns a chebfun2 representing the text in
%  STRING. The full English alphabet is supported.
%
%  Example:
%   f = scribble('Happy birthday LNT!');
%   contour(f,.2:.1:.85), axis equal
%
% See also scribble.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )
    s = 'HAPPY BIRTHDAY LNT!';
end

s = upper(s);

space = zeros(9,9);

A = space;
A(1,[4,5,6]) = 1;
A(2,[3,4,5,6,7]) = 1;
A(3,[2,3,4,6,7,8]) = 1;
A(4,[1:3,7:9]) = 1;
A(5:end,[1:3,7:9]) = 1;
A([5,6],:) = 1;

B = space;
B([1,2,4,5,8,9],:) = 1;
B(:,1:2) = 1;
B(:,8:9) = 1;
B(81) = 0;
B(73) = 0;
B(77) = 0;

C = space;
C([1,9],4:8) = 1;
C([2,8],[3:5 7:9]) = 1;
C([3,7],[2:4 8:9]) = 1;
C([4,5,6],1:3) = 1;

D = space;
D(:,[1,2,8,9]) = 1;
D([1,2,8,9],:) = 1;
D([73, 81]) = 0;

E = space;
E(:,[1,2]) = 1;
E([1,2,4,5,6,8,9],2:9) = 1;
E(4:7,8:9) = 0;

F = E;
F(8:9,3:9) = 0;
F(6,3:9) = 0;
F(4:5,7) = 0;

G = space;
G([1,9],3:8) = 1;
G([2,8],2:9) = 1;
G([3,7],[1:3 7:9]) = 1;
G([4,5,6],1:2) = 1;
G([5,6],6:9)=1;

H = space;
H(:,[1,2,8,9]) = 1;
H([4,5],:) = 1;

I = space;
I(:,4:6) = 1;
I([1,2,8,9],1:9) = 1;

J = space;
J([1,2],1:end) = 1;
J(1:end-2,[4,5,6]) = 1;
J(7,1:2) = 1;
J(8,1:6) = 1;
J(9,2:5) = 1;

K = space;
K(:,[1,2]) = 1;
K(4:6,3) = 1;
K(32:10:end) = 1;
K(33:10:63) = 1;
K(31:10:end) = 1;
K(32:8:70) = 1;
K(31:8:70) = 1;
K(33:8:end) = 1;

L = space;
L(:,[2,3]) = 1;
L([8,9],2:9) = 1;

M = space;
M(:,[1,2,8,9]) = 1;
M(20:10:40) = 1;
M(21:10:41) = 1;
M(41:8:71) = 1;
M(40:8:71) = 1;

N = space;
N(:,[1,2,8,9]) = 1;
N(10:10:end) = 1;
N(11:10:end) = 1;
N(12:10:end) = 1;

O = space;
O([1,2,8,9],:) = 1;
O(:,[1,2,8,9]) = 1;
O([1,9,73,81]) = 0;

P = space;
P(:,[1,2]) = 1;
P([1,2,5,6],:) = 1;
P([3,4],8:9) = 1;
P(1,9) = 0;
P(6,9) = 0;

Q = O;
Q(end) = 1;
Q(60:61) = 1;
Q(51:52) = 1;
Q(79) = 0;
Q(63) = 0;

R = P;
R(32:10:end) = 1;
R(41:10:end) = 1;
R(23:10:end) = 1;
R(1,9) = 0;

S = space;
S([1,2,4,5,8,9],:) = 1;
S(1:5,1:2) = 1;
S(5:9,8:9) = 1;
S(1,1) = 0;
S(end) = 0;
S(4,9) = 0;
S(5,1) = 0;

T = space;
T(1:2,:) = 1;
T(:,4:6) = 1;

U = O;
U(1:2,:) = 0;
U(1:2,[1,2,8,9]) = 1;

W = space;
W(1:12:30) = 1;
W(2:12:30) = 1;
W(3:12:30) = 1;
W(11:12:30) = 1;
W(12:12:30) = 1;
W(10:12:30) = 1;
W(8:9,4) = 1;
W(7:8,5) = 1;
W(7,4:6) = 1;
W(6,5) = 1;
W = double(W | fliplr(W));

V = space;
V(1:11:50) = 1;
V(2:11:40) = 1;
V(3:11:40) = 1;
V(7:8,5) = 1;
V(1:2,2) = 1;
V(4,3) = 1;
V(6,4) = 1;
V(:,5:9) = fliplr(V(:,1:5));

X = space;
X(1:10:end) = 1;
X(2:10:end) = 1;
X(10:10:end) = 1;
X = double(X | fliplr(X));

Y = space;
Y(1:10:40) = 1;
Y(2:10:40) = 1;
Y(10:10:40) = 1;
Y(5:9,4:5) = 1;
Y = double(Y | fliplr(Y));

Z = space;
Z([1:2,8:9],:) = 1;
Z(9:8:end) = 1;
Z(8:8:end) = 1;
Z(10:8:end) = 1;

%Numerals
zero = O;
zero(5,5) = 1;

one = I;
one(1:2,[1:3 7:9]) = 0;
one(2,3) = 1;
one(3,2:3) = 1;
one(4,1:2) = 1;

two = space;
two(1,2:8) = 1;
two(2,:) = 1;
two(3,[1:2 7:9]) = 1;
two(4,6:8) = 1;
two(5,5:7) = 1;
two(6,4:6) = 1;
two(7,3:5) = 1;
two(8:9,:) = 1;

three = fliplr(E);
three([4 6],3:6) = 0;
three([1 5 9],9) = 0;
three([3 7],7) = 1;
three(5,2) = 1;

four = space;
four(:,6:8) = 1;
four(6:7,:) = 1;
four(2,5) = 1;
four(3,4:5) = 1;
four(4,3:4) = 1;
four(5,2:3) = 1;

five = S;
five([1 5],1) = 1;
five(7,1:2) = 1;
five(9,1) = 0;

six = S;
six(5:8,1:2) = 1;
six(9,1) = 0;
six(1,9) = 0;

seven = two;
seven(1,[1 9]) = 1;
seven(8,[1 5:9]) = 0;
seven(9,4:9) = 0;

eight = B;
eight([1 5 9]) = 0;

nine = S;
nine(3:5,8:9) = 1;
nine([1 9],[1 9]) = 0;
nine(7,1:2) = 1;

% Punctuation
exclam = space;
exclam([1:6 8:9],4:6) = 1;

% Padding:
b = zeros(9,2);
Str = zeros(9,50);
pad = [];

for jj = 1:length(s)
    Str = [Str ; zeros(1, size(Str, 2))]; %#ok<AGROW>
    b = [b ; 0, 0];                 %#ok<AGROW>
    pad = [pad ; zeros(1, 9)];      %#ok<AGROW>
    
    switch s(jj)
        case {' '}
            Str = [ [Str, b], [pad ; space]];
        case {'!'}
            Str = [ [Str, b], [pad ; exclam]];
        case {'0'}
            Str = [ [Str, b], [pad ; zero]];
        case {'1'}
            Str = [ [Str, b], [pad ; one]];
        case {'2'}
            Str = [ [Str, b], [pad ; two]];
        case {'3'}
            Str = [ [Str, b], [pad ; three]];
        case {'4'}
            Str = [ [Str, b], [pad ; four]];
        case {'5'}
            Str = [ [Str, b], [pad ; five]];
        case {'6'}
            Str = [ [Str, b], [pad ; six]];
        case {'7'}
            Str = [ [Str, b], [pad ; seven]];
        case {'8'}
            Str = [ [Str, b], [pad ; eight]];
        case {'9'}
            Str = [ [Str, b], [pad ; nine]];
        otherwise
            Str = [ [Str, b], [pad ; eval(s(jj))]];
    end
end

% Padding:
Str = [Str, zeros(size(Str, 1), 50)];
Str = [zeros(20, size(Str, 2)); Str ; zeros(20, size(Str, 2))];

if (nargin == 2)
    f = chebfun2( flipud(Str), rk );
else
    f = chebfun2( flipud(Str) );
end

% if ( nargout == 0 )
%    % Plot: 
%     contour(f,.3:.1:.85,'numpts',1000), axis off
%     rk = length(f);
%     t = sprintf('Rank = %u',rk);
%     title(t,'fontsize',16)
% else
%     varargout = {f};
% end
x = linspace(-1,1,200); 
[xx,yy] = meshgrid(x); 
[C, D, R] = cdr( f ); 
C = C(x',:); 
R = R(x',:); 
if ( nargout == 0 )
    r = rank(f);
    for rk = 1:r
        use = rk; 
        if ( rk == 11 ) 
            use = 10; 
        end
        if ( rk == 16 ) 
            use = 15; 
        end
        g = C(:,1:use)*D(use,use)*R(:,1:use)'; %chebfun2(flipud(Str), rk);
        contour(xx,yy,g,.35:.05:.85), axis off
%        set(gcf, 'color', 'w')
        %rk = length(f);
        t = sprintf('Rank = %u',rk);
        title(t,'fontsize',16), shg
%         im = frame2im(getframe());
%         [imind, cm] = rgb2ind(im, 16);
%         if ( rk == 1 )
%             imwrite(imind, cm, 'Scribble2.gif', 'gif', ...
%                 'Loopcount', inf, 'DelayTime', 1e-5);
%         else
%             imwrite(imind, cm, 'Scribble2.gif', 'gif', ...
%                 'WriteMode', 'append', 'DelayTime', 1e-5);
%         end
        shg
        drawnow
        
    end
else
 	varargout = {f}; 
end




