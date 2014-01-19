%% Event handling in Chebfun2v 
% Alex Townsend, 9th June 2013

%% 
% (Chebfun2 example veccalc/EventHandling.m)
% [Tags: #Bouncing, #autonomous, #Chebfun2]

function EventHandling
%%
LW = 'linewidth'; lw = 2; 
FS = 'fontsize'; fs = 16;
MS = 'markersize'; ms = 20;

%%
% A few days ago Jonathan Black from the OCIAM group in Oxford University
% was interested in using the chebfun2v command ODE45 for solving an
% autonomous system with event handling. Unfortunately, Chebfun2 did not
% provide support for event handling. This Example details two small
% features added to the ODE45 command in Chebfun2: 

%%
% 1. Ability to solve autonomous systems with Chebfun2v objects with three
% components; and 

%%
% 2. Event handling. That is, chebfun2's ODE45 can now be used with the
% standard Matlab's ODE45 event handling features.

%%
% For a different example with Chebfun2's ODE45 command see [2] and the
% Chebfun2 guide.

%%
% In principle Chebfun ODE45 might have events too, though at present it 
% doesn't. This may change in the future. 

%%
% We compute the trajectory of a ball bouncing down stairs using the
% new features of Chebfun2's ODE45 command. For a movie of a ball bouncing 
% see [1]. Let's make it interesting by representing the stairs by a 
% piecewise chebfun: 

steps = @(t) 1.5*(0<=t & t<2) + (2<=t & t<4) + .5*(4<=t & t<6) ... 
      + .25*(6<=t & t<8) + .1*(8<=t & t<10) + .5*(10<=t & t<12);
stairs = chebfun(steps, [0,30], 'splitting', 'on');
plot(stairs,LW,lw), title('Stairs',FS,fs), axis([0 25 0 2])

%%
% Also let's suppose the ball is initially at 2 metres above the floor, 
% $h(0)=2$, and is thrown horizontally at a speed of one metre per second,

u0 = [2 0 0];  % u0 = [h(0) h'(0) x(0)]

%% 
% Further let's suppose the ball then bounces down the stairs under 
% Martian gravity (gravity constant of 1) with a small Martian air resistence. 
% We have the following equations to model the height of the ball $h(t)$:
% 
%  $$h''(t) = -1-.01h',\quad h(0) = 2,\quad h'(0) = 0,$$ 
%
% and for the horizontal distance from the origin we have 
% 
%  $$x'(t) = 1,\quad x(0) = 0.$$
% 
% This can be rewritten as an autonomous system by introducing the dummy
% variable $v(t) = h'(t)$, and we can represent the system in the following
% chebfun2v object (see [2]):  

t0 = 0; Tend = 50;  
F = chebfun2v(@(h,hp)hp, @(h,hp)-1-.01*hp, @(h,hp)1+0*h,[0 30 0 2]);

%% 
% Once the ball makes contact with the ground we suppose that the ball
% loses 15 percent of its energy.  We need
% event handling in the Chebfun2's ODE45 so when the ball
% bounces we can reverse the vertical velocity of the ball. 
% We terminate the ODE45 solver when the maximum height of a bounce
% is below 5cm.  

%%
% Note we are using event handling to simulate a slabwise continuous 
% autonomous system in 3D, which cannot be represented by a chebfun2v on
% its own because a chebfun2v is a smooth vector field. 

%%
% Here is the path taken by the bouncing ball:

% Event when the ball hits the ground. 
options = odeset('RelTol',100*eps, 'events', @p1Event1);
h = []; v = []; x = [];   % store solution
while ( 1 )  
    sol = ode45(F, [t0, Tend], u0, options); 
    t0 = sol.xe;
    h = [h; sol.y(:,1)]; x = [x; sol.y(:,3)];
    if max(sol.y(:,1)) < 5e-3, break, end
    u0 = (sol.ye).*([1;-.85;1]); %  reverse the velocity
end
stairs = stairs{0,max(x)};
plot(stairs,LW,2), hold on, plot(x,h,LW,lw)
grid on, xlabel('time')
ylabel('h(t)'), axis([stairs.ends([1 end]) 0 2])
title('Bouncing ball',FS,fs)

%% 
% The path of the bouncing ball is now represented by a chebfun and we can 
% therefore work out the maximum height of each bounce with Chebfun's roots
% command. Since we do not want to count the bounces we supply the roots
% command with the flag 'nojump'.

d = roots(diff(h) ,'nojump'); 
plot(x(d),h(d),'.r',MS,ms)

end

function [value, isTerminal, direction] = p1Event1(t,y)
% Event handling function. Stop Chebfun2's ODE45 when the ball bounces. 
steps = @(x) 1.5*(0<=x & x<2) + (2<=x & x<4) +... 
    .5*(4<=x & x<6) + .25*(6<=x & x<8) + .1*(8<=x & x<10) + .5*(10<=x & x<12);
value = y(1)-steps(y(3)); 
isTerminal = 1; 
direction = -1; 
end

%%  
% References
%%
% [1] F. Di Tommaso, A bouncing ball, Chebfun Example, 
% http://www2.maths.ox.ac.uk/chebfun/examples/ode/html/BouncingBall.shtml,
% February 2013.
%%
% [2] A. Townsend, Phase portraits and trajectories, Chebfun2 Example,
% http://www2.maths.ox.ac.uk/chebfun/examples/veccalc/html/AutonomousSystems.shtml,
% March 2013. 