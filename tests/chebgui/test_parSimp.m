function pass = test_parSimp(varargin)

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

S = loadtests;

pass = zeros(numel(S), 1);
for k = 1:numel(S)
    str = stringParser.parSimp(S{k}{1});
    pass(k) = strcmp(str,S{k}{2});
end


function S = loadtests
% BVPs
S{1} = {'((0.02.*diff(u,2)+diff(u))+u)',
    '0.02.*diff(u,2)+diff(u)+u'};
S{2} = {'((0.01.*diff(u,2)-x.*u)-1)',
    '0.01.*diff(u,2)-x.*u-1'};
S{3} = {'((x.^(2).*diff(u,2)+x.*diff(u))+(x.^(2)-3.^(2)).*u)',
    'x.^2.*diff(u,2)+x.*diff(u)+(x.^2-3.^2).*u'};
S{4} = {'(((0.01.*diff(u,2)+2.*(1-x.^(2)).*u)+u.^(2))-1)',
    '0.01.*diff(u,2)+2.*(1-x.^2).*u+u.^2-1'};
S{5} = {'((diff(u,2)+(1.2+sign((10-abs(x)))).*u)-1)',
    'diff(u,2)+(1.2+sign((10-abs(x)))).*u-1'};
S{6} = {'((diff(u,2)+u)-u.^(2))',
    'diff(u,2)+u-u.^2'};
S{7} = {'(diff(u,4)-(diff(u).*diff(u,2)-u.*diff(u,3)))',
    'diff(u,4)-diff(u).*diff(u,2)+u.*diff(u,3)'};
S{8} = {'(diff(u,2)+.87.*exp(u))',
    'diff(u,2)+.87.*exp(u)'};
S{9} = {'((.0005.*diff(u,2)+x.*(x.^(2)-0.5).*diff(u))+3.*(x.^(2)-0.5).*u)',
    '.0005.*diff(u,2)+x.*(x.^2-0.5).*diff(u)+3.*(x.^2-0.5).*u'};
S{10} = {'((.01.*diff(u,2)+u.*diff(u))-u)',
    '.01.*diff(u,2)+u.*diff(u)-u'};
S{11} = {'(diff(u,2)+sin(u))',
    'diff(u,2)+sin(u)'};
S{12} = {'((0.05.*diff(u,2)+diff(u).^(2))-1)',
    '0.05.*diff(u,2)+diff(u).^2-1'};
S{13} = {'(diff(u,2)-8.*sinh(8.*u))',
        'diff(u,2)-8.*sinh(8.*u)'};
S{14} = {'(diff(u,2)-(u-1).*(1+diff(u).^(2)).^(1.5))',
    'diff(u,2)-(u-1).*(1+diff(u).^2).^1.5'};
S{15} = {'((diff(u,2)-x.*sin(u))-1)',
    'diff(u,2)-x.*sin(u)-1'};
S{16} = {'((diff(u,2)-(1-u.^(2)).*diff(u))+u)',
    'diff(u,2)-(1-u.^2).*diff(u)+u'};
S{17} = {'(diff(u,2)-sin(v))',
    'diff(u,2)-sin(v)'};
S{18} = {'(diff(u,2)-sin(v))',
    'diff(u,2)-sin(v)'};
S{19} = {'(cos(u)+diff(v,2))',
    'cos(u)+diff(v,2)'};

% PDEs
S{20} = {'+(0.1.*diff(u,2)+diff(u))',
    '0.1.*diff(u,2)+diff(u)'};
S{21} = {'+((.01.*diff(u,2)+u)-u.^(3))',
    '.01.*diff(u,2)+u-u.^3'};
S{22} = {'+(-diff(u.^(2))+.02.*diff(u,2))',
    '-diff(u.^2)+.02.*diff(u,2)'};
S{23} = {'+(-.003.*diff(u,4)+diff((u.^(3)-u),2))',
    '-.003.*diff(u,4)+diff(u.^3-u,2)'};
S{24} = {'+0.1.*diff(u,2)',
    '0.1.*diff(u,2)'};
S{25} = {'+(.02.*diff(u,2)+cumsum(u).*sum(u))',
    '.02.*diff(u,2)+cumsum(u).*sum(u)'};
S{26} = {'+((u.*diff(u)-diff(u,2))-0.006.*diff(u,4))',
    'u.*diff(u)-diff(u,2)-0.006.*diff(u,4)'};

% Systems
S{27} = {'+((-u+(x+1).*v)+0.1.*diff(u,2))',
    '-u+(x+1).*v+0.1.*diff(u,2)'};
S{28} = {'+((u-(x+1).*v)+0.2.*diff(v,2))',
    'u-(x+1).*v+0.2.*diff(v,2)'};
S{29} = {'+(0.1.*diff(u,2)-100.*u.*v)',
    '0.1.*diff(u,2)-100.*u.*v'};
S{30} = {'+(.2.*diff(v,2)-100.*u.*v)',
    '.2.*diff(v,2)-100.*u.*v'};
S{31} = {'+(0.001.*diff(w,2)+200.*u.*v)',
    '0.001.*diff(w,2)+200.*u.*v'};
S{32} = {'+(diff(u,2)-v)',
    'diff(u,2)-v'};
S{33} = {'(diff(v,2)-u)',
    'diff(v,2)-u'};
S{34} = {'+(diff(u,2)+exp(-100.*x.^(2)).*sin(pi.*t))',
    'diff(u,2)+exp(-100.*x.^2).*sin(pi.*t)'};

% EIGs
S{35} = {'(diff(u,2)+diff(u))',
    'diff(u,2)+diff(u)'};
S{36} = {'(-diff(u,2)+1i.*x.^(2).*u)',
    '-diff(u,2)+1i.*x.^2.*u'};
S{37} = {'(-.1.*diff(u,2)+4.*(sign((x+1))-sign((x-1))).*u)',
    '-.1.*diff(u,2)+4.*(sign((x+1))-sign((x-1))).*u'};
S{38} = {'(-diff(u,2)+x.^(2).*u)',
    '-diff(u,2)+x.^2.*u'};
S{39} = {'(-diff(u,2)+5.*cos(2.*x).*u)',
    '-diff(u,2)+5.*cos(2.*x).*u'};
S{40} = {'((diff(u,2)+u.*x)+v)',
    'diff(u,2)+u.*x+v'};
S{41} = {'(diff(v,2)+sin(x).*u)',
    'diff(v,2)+sin(x).*u'};

% MISC

% Check exponential notation
S{42} = {'(diff(v,2)-1e-2*u)',
    'diff(v,2)-1e-2*u'};
S{43} = {'(diff(v,2)-1e+2*u)',
    'diff(v,2)-1e+2*u'};
S{44} = {'(diff(v,2)+1e-2*u)',
    'diff(v,2)+1e-2*u'};
S{45} = {'(diff(v,2)+1e+2*u)',
    'diff(v,2)+1e+2*u'};