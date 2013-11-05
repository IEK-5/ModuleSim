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
Rp=1e8;				% parallel resistance of the junction				[Ohm cm2]
Js=1e-8;			% Saturation current						[A / cm2]
Jph=0.03;			% Photo current							[A / cm2]
nid=1.5;			% ideality factor						[-]
T=300;				% Temperature							[K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization:
% we (must!) first define the diode 
VJ = GenJV(Rs, Rp,  Js, Jph, nid, T, -0.1, 1.0, 100);	% This produces a diode JV characteristic including series and parallel resistance
SetDiode(VJ)						% this command sets the diode JV characteristic to the system

%%%%%%%%%%%%%% insert shunts (at least one shunt)
% We define a simple linear JV characteristic:
ohm=10; 					% The total resistance over shunted area
r=0.02; 					% The radius of the shunt radius
rr=pi*r^2*ohm; 					% Compute the area resistance in Ohm cm^2
VJ=[-1,-1/rr;1,1/rr]; 				% the resulting JV characteristic
AddPointSource(VJ,(x1+x2)/2,(y1+y2)/2, r) 	% add the point source to the system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute an VI curve
N=20;
Vstart=0;
Vstop=0.6;
res=[];
for i=0:N
	V=Vstart+i*(Vstop-Vstart)/N;
	[I,Iext]=SetVoltage(V, 1e-3); 	% This command solves the operating point
	res=[res;V,sum(Iext)];		% The externbal current is an array of the external currents induced by each part-solution
					% it's sum is thus the total current.
	printf("%e %e\n",V,sum(Iext));
endfor
SaveXY("IV.dat",res);
