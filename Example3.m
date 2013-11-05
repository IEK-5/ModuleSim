#!/usr/bin/octave -qf
source ThinFilmSim.m;
warning("on", "backtrace")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters:
x1=0;				% x coordinate for the lower left corner 			[cm]
x2=0.4;				% x coordinate for the upper right corner 			[cm]
y1=0;				% y coordinate for the lower left corner 			[cm]
y2=20;				% y coordinate for the upper right corner 			[cm]
Rf=18;				% Sheet resistance of the front contact 			[Ohm]

Rs=1e-3;			% series resistance of the junction				[Ohm cm2]
Rp=1e3;				% parallel resistance of the junction				[Ohm cm2]
Js=1e-8;			% Saturation current						[A / cm2]
Jph=0.0;			% Photo current							[A / cm2]
nid=1.5;			% ideality factor						[-]
T=300;				% Temperature							[K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization:

% we (must!) first define the diode 
VJ = GenJV(Rs, Rp,  Js, Jph, nid, T, -0.1, 1.0, 100);	% This produces a diode JV characteristic including series and parallel resistance
SetDiode(VJ)						% this command sets the diode JV characteristic to the system

% We define a simple linear JV characteristic:
ohm=1; 						% The total resistance over shunted area
r=0.02; 					% The radius of the shunt radius
rr=pi*r^2*ohm; 					% Compute the area resistance in Ohm cm^2
VJ=[-1,-1/rr;1,1/rr]; 				% the resulting JV characteristic
AddPointSource(VJ,(x1+x2)/2,(y1+y2)/2, r) 	% add the point source to the system

% insert an illuminated spot
VJ=GenJV(Rs, Rp, Js, 0.03, nid, T, -0.1, 1.0, 100);
r=0.0010;
AddPointSource(VJ,0.0,0.0, r);% for now the position is not important

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate an LBIC experiment
% generate a mesh, we will position the illuminated spot one by one on all the coordinates of this mesh
% The next block produces a list of coordinates (c), every row is one coordinate where the first column contains 
% the x values and the second the y values

Nx=20;		% resolution = 0.4/Nx
Ny=500;		% resolution = 20/Ny

x=x1+0.5*(x2-x1)/Nx:(x2-x1)/Nx:x1+(Nx-0.5)*(x2-x1)/Nx;
y=9+0.5*(y2-y1)/Ny:(y2-y1)/Ny:11-0.5*(y2-y1)/Ny;

NN=length(x)*length(y);
[xx,yy]=meshgrid(x,y);
x=reshape (xx,NN, 1);
y=reshape (yy,NN, 1);
c=[x,y];


% now we can start: first initialize our result varable
res=[];
% it is going to take a while so I will present you with a progress bar:
pco=-1;	% old percentage complete (needs to be initialized negative)
Nc=length(c(:,1));
for i=1:Nc
	SetPointSource(2,VJ,c(i,1),c(i,2),r)		% change the position of the illuminated spot (source with index 2) to the next coordinate
	[I,Iext]=SetVoltage(0.0, 1e-6); 		% Solve the operating point (short circuit)
	res=[res;c(i,:),sum(Iext), Iext'];		% store the results
	pcn=100*i/Nc;					% new percentage complete
	pco=ProgressBar(pcn, pco, 40, 4); 		% draw the progress bar
endfor
SaveXY("lbic.dat",res);
