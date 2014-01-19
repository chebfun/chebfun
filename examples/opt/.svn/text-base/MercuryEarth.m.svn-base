%% Mercury-Earth minimum separation
% Tonatiuh Sanchez-Vizuet and Matthew Moye, 19th June 2012

%%
% (Chebfun example opt/MercuryEarth.m)
% [Tags: #optimization, #astrophysics]

%%
% The Earth and Mercury on their elliptical orbits around the sun can be
% shown to have minimized separation at a time t. 
%
% The domain will be the time duration of evaluation in days.
domain = [0 1000];
t = chebfun('t', domain);

%%
% The parametrized equations for the orbits are given by [1]:
y_m = 56.6741*sin(2*pi*t/87.97); 
x_m = -11.9084+57.9117*cos(2*pi*t/87.97);
y_e = 149.5832*sin(2*pi*t/365.25);
x_e = -2.4987 + 149.6041*cos(2*pi*t/365.25); 
%%
% Chebfun is excellent in taking a function like the distance equation over
% a given domain.
f = sqrt((y_m-y_e).^2 + (x_m-x_e).^2);

%%
% Minival is the minimum distance and minpos is the time of occurrence. 
[minival,minpos] = min(f);

plot(t,f)
xlabel('Time (days)')
hold on, plot(minpos,minival, '.r', 'markersize', 20)

%%
% First, the parametrized curves are plotted.
figure
 hold off
 plot(x_m, y_m), hold on
 plot(x_e, y_e)
 
%% 
% Orbital depiction of the planets' positions at minpos.
 plot(x_m(minpos),y_m(minpos),'.r', 'markersize', 20)
 plot(x_e(minpos),y_e(minpos),'.r', 'markersize', 20)
 title('Mercury and Earth Orbits')
 
%%
% References:
%
% [1] Charles F. Van Loan, Introduction to Scientific Computing, 1997, p.
% 274.
