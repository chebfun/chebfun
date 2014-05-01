function x = fourierpts(n)

x = linspace(-pi, pi, n+1).';
x = x./pi;
x(end) = [];

end