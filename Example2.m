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
% We define a simple linesr JV characteristic:
ohm=10; 					% The total resistance over shunted area
r=0.02; 					% The radius of the shunt radius
rr=pi*r^2*ohm; 					% Compute the area resistance in Ohm cm^2
VJ=[-1,-1/rr;1,1/rr]; 				% the resulting JV characteristic
AddPointSource(VJ,(x1+x2)/2,(y1+y2)/2, r) 	% add the point source to the system

% I add another two defects:
ohm=5;
r=0.02;
rr=pi*r^2*ohm;
VJ=[-1,-1/rr;1,1/rr];
AddPointSource(VJ,(x1+x2)/2,(y1+y2)/4, r)	% note the y coordinate is changed w.r.t. the first defect
 
ohm=20;
r=0.02;
rr=pi*r^2*ohm;
VJ=[-1,-1/rr;1,1/rr];
AddPointSource(VJ,(x1+x2)/2,3*(y1+y2)/4, r)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define a single operating point:
% here we inject a current of -200 mA
[I, Ve, Iext]=SetCurrent(-0.2, 1e-6)

% alternatively one could set an external voltage:
% [I,Iext]=SetVoltage(0.4, 1e-6)

% or set the maximum powerpoint as the operating point
% [I,V, Iext]=FindMPP(1e-4)

% Now let us get some data out:
% we start by defining a mesh:
Nx=40;
Ny=200;

x=x1+0.5*(x2-x1)/Nx:(x2-x1)/Nx:x1+(Nx-0.5)*(x2-x1)/Nx;
y=y1+0.5*(y2-y1)/Ny:(y2-y1)/Ny:y2-0.5*(y2-y1)/Ny;
NN=length(x)*length(y);
[xx,yy]=meshgrid(x,y);
x=reshape (xx,NN, 1);
y=reshape (yy,NN, 1);
c=[x,y];

% the variable c now consists of two columns, one for x coordinates and one for y coordinates
% The function VEJP computes the potential (V) the electric field (E) the current densities (J) and the power densities (P)
% as input is takes a an array of coordinates and the current vector with the current for each part solution (as determined
% above).
[V,Ef,Jf, Jj, Pf,Pj]=VEJP(c,I);
SaveXY("V.dat",V);
SaveXY("Ef.dat",Ef);
SaveXY("Jf.dat",Jf);
SaveXY("Jj.dat",Jj);
SaveXY("Pf.dat",Pf);
SaveXY("Pj.dat",Pj);
% Note that the power density within the point defects is not computed in VEJP (the values within the defect are indeed still the
% power densities dissipated in the junction, the point defect is in parallel to the junction.
% However, we already have the current flowing through each part solution in I, which includes all the point defects. Furthermore we
% have all the JV characteristics of each solution. Using this we can easily compute the total power dissipated in the diode and in each
% defect:
Vp=SolV(I);
Pp=Vp.*I
