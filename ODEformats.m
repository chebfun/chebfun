set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
set(0, 'defaultfigureposition', [380 320 540 200]);
set(0, 'defaultaxeslinewidth',  0.7);
set(0, 'defaultaxesfontsize',   7);
set(0, 'defaultlinelinewidth',  .9);
set(0, 'defaultpatchlinewidth', .9);
set(0, 'defaultlinemarkersize', 15); 
set(0, 'defaultaxesfontweight', 'normal'); 
set(0, 'defaulttextinterpreter', 'latex'); 
format compact
format short
chebfunpref.setDefaults('factory');
FS = 'fontsize'; LW = 'linewidth'; MS = 'markersize'; CO = 'color';
IN = 'interpret'; LT = 'latex';
XT = 'xtick'; YT = 'ytick';
XTL = 'xticklabel'; YTL = 'yticklabel';
LO = 'location'; NE = 'northeast'; NO = 'north';
HA = 'HorizontalAlignment'; CT = 'center'; RT = 'right';
FN = 'fontname'; YS = 'ystretch'; LS = 'linestyle';
purple = [.8 0 1]; green = [.466 .674 0]; %green = [0 .7 0];
blue = [0 .447 .741];
%ivp = [.15 .8 0]; ivpnl = [0 .35 0];
ivp = [.466 .674 0]; ivpnl = [.23 .34 0];
%bvp = [0 0 1]; bvpnl = [0 0 .5];
bvp = [0 .447 .741]; bvpnl = [0 .23 .37];
ibvp = [.85 0 .8]; ibvp0 = [.5 0 .4];
orange = [1 .5 0];
ibvp = orange; ibvp0 = .6*ibvp;
