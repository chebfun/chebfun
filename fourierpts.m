function x = fourierpts(n, dom)


x = linspace(-pi, pi, n+1).';
x = x./pi;
x(end) = [];

if ( nargin == 1 )
    dom = [-1, 1]; 
end

% map to the right domain: 
x = diff(dom(1:2))*x/2 + mean(dom(1:2)); 
    
end