#!/usr/bin/octave -qf
% Model for a cell stripe with shunts using the sourcefield program to compute the influence of the shunts
% The model includes the following parameters and assumptions
% contacts are equipotential
% a perfect back contact is assumed
% the front contact sheet resistance
% the non-linear junction JV is linearized around the operating point (one operating point for the entire domain)
% non-linear shunts, are also linearized around the operating point

if (! exist ("SourceField", "file"))
	fprintf(stderr,"Error: The SourceField executable is not found\n");
	fprintf(stderr,"       Please place a copy of the SourceField program in:\n       %s \n",  pwd ());
	exit();
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters:
global x1=0;				% x coordinate for the lower left corner 			[cm]
global x2=0.4;				% x coordinate for the upper right corner 			[cm]
global y1=0;				% y coordinate for the lower left corner 			[cm]
global y2=20;				% y coordinate for the upper right corner 			[cm]
global Rf=18;				% Sheet resistance of the front contact 			[Ohm]

global VJR={};				% variable to be filled by script, contains the VJ tables + the resistance
global lin_sys=[];			% linearisation of the system
global Nbisect=3;			% every Nbisect iterations one more robust iteration to force convergence
global Nnumint=500;			% number op points to compute the integral of the external cuurent induced by a local source
global SFErr=1e-4;			% the error for the sourcefield calculation
global k=8.617e-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Linearize and generate JV curves
% Lambert-W function
function [y,n]=W(x);
	% first estimate
	y=log(1+x);
 	y=y.*(1.-(log(1+y))./(2.+y));
	N=50; % maximum number of iterations
	e=0.5e-15; % acc to about machine precision
	n=0; % number of iterations
	E=abs(x.-y.*exp(y))./(abs(x)+e);
	while ((n<N) && (max(E)>e))
		n++;
		y=(y-(y.*exp(y)-x)./(exp(y).*(y+1)-(y+2).*(y.*exp(y)-x)./(2.*y+2)));
		E=abs(x.-y.*exp(y))./(abs(x)+e);
	endwhile
	[m,im]=max(E);
	if (m>e)
		fprintf(stderr, "Warning: Iteration limit reached in LambertW evaluation\n");
		[m,im]=max(E);
		fprintf(stderr, "Maximum relative error estimated to be %e for x=%e\n", m, x(im));
	endif
endfunction

function J=Diode(V, Rs, Rp, Js, Jph, nid, T)
	global k;
	nkT=nid*k*T;	
	J=(nkT/Rs)*W((Rs*Rp*Js/((Rs+Rp)*nkT))*exp(((Rp/(Rs+Rp))*V+(Jph+Js)*Rp*Rs/(Rp+Rs))/nkT))-((Jph+Js)*Rp-V)/((Rs+Rp));
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Linearize and generate JV curves
% Generates a simple diode equation JV curve
function JV = GenJV(Rs, Rp,  Js, Jph, nid, T, V1, V2, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	Rs: 	Series resistance 	[Ohm cm2]
%	Rp: 	Parallel resistance 	[Ohm cm2]
%	nd: 	Diode ideality factor 	[-]
%	Js: 	Saturation Current 	[A/cm2]
%	nd: 	Diode ideality factor 	[-]
%	T: 	Temperature 		[K]
%	Jph: 	Photocurrent 		[A/cm2]
%	V1: 	Start voltage 		[V]
%	V2: 	End voltage 		[V]
%	N: 	Number of steps 	[-]
% OUTPUT:
%	JV:	JV curve in two columns	[V]:[A/cm2]
	JV=[];
	Vstep=(V2-V1)/(N-1);
	V=V1:Vstep:V2;
	J=J=Diode(V, Rs, Rp, Js, Jph, nid, T);
	JV=[V',J']; % all data is there, now we have to reorder the array to something like this:
endfunction 

% differentiates a single column
function dX=Diff(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	X: 	column vector 		[?]
% OUTPUT:
%	dX:	dX			[d?]
	X=X(:,1);
	dX=X;
	dX(2:end-1)=(X(3:end).-X(1:end-2))./2;
	dX(1)=X(2)-X(1);
	dX(end)=X(end)-X(end-1);
endfunction


% Next to functions linearize the junction as a function of voltage and current
function [Rj,Vc]=LinVI_V(V,VJR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	V: 	Voltage 		[V]
%	VJR:	columnar data		[V]:[A/cm2]:[Ohm cm2]
% OUTPUT:
%	Rj:	Junction resistance	[Ohm cm2]
%	Vc:	Junction offset voltage	[V]
	if ((V>VJR(1,1)) && (V<VJR(end,1)))
		J=interp1(VJR(:,1),VJR(:,2),V);
		Rj=interp1(VJR(:,1),VJR(:,3),V);
		Vc=V-J*Rj;
	else
		if(V<=VJR(1,1))
			J=VJR(1,2);
			Rj=VJR(1,3);
			Vc=VJR(1,1)-J*Rj;
		else
			J=VJR(end,2);
			Rj=VJR(end,3);
			Vc=VJR(end,1)-J*Rj;
		endif
	endif
endfunction

function [Rj,Vc]=LinVI_I(J,VJR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	J: 	current density		[A/cm2]
%	VJR:	columnar data		[V]:[A/cm2]:[Ohm cm2]
% OUTPUT:
%	Rj:	Junction resistance	[Ohm cm2]
%	Vc:	Junction offset voltage	[V]
	if ((J>VJR(1,2)) && (J<VJR(end,2)))
		V=interp1(VJR(:,2),VJR(:,1),J);
		Rj=interp1(VJR(:,2),VJR(:,3),J);
		Vc=V-J*Rj;
	else
		if(J<=VJR(1,2))
			V=VJR(1,1);
			Rj=VJR(1,3);
			Vc=V-VJR(1,2)*Rj;
		else
			V=VJR(end,1);
			Rj=VJR(end,3);
			Vc=V-VJR(end,2)*Rj;
		endif
	endif
endfunction			
			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute superimposed potantial solutions
% Write a mesh for sourcefield, c is an array with coordinates [x1,y1;x2,y2;...;xN,yN]
function WriteMesh(c,fn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
%	fn:	file name		[-]
% OUTPUT:
%	Columnar data written tab separated to plain ascii file
	f=fopen(fn, "w");
	for i=1:length(c(:,1))
		fprintf(f,"%e\t%e\n", c(i,1), c(i,2));
	endfor	
	fclose(f);
endfunction

% Computes the potential acording to the formula for a cell without shunts (i.e. 1 1d model only depends on x value)
% The potential is normalized to an injected current of 1A
function [V, Ef]=BasePotential(c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
% OUTPUT:
%	V:	Voltage @ c		[V]
%	Ef:	Electric Field in x&y	[V/cm]:[V/cm]
%		direction
	global Rf;
	global lin_sys;
	global y1;
	global y2;
	global x1;
	global x2;
	Rj=lin_sys(1,4);
	
	rj=Rj;
	rf=Rf;
	lambda=sqrt(rf/rj);
	Jo=1/(y2-y1);
	Vo=rf*Jo*cosh(lambda*(x2-x1))/(lambda*sinh(lambda*(x2-x1)));
	
	% compute potential
	V=Vo*cosh(c(:,1).*lambda)-Jo*rf*sinh(c(:,1).*lambda)/lambda;
	
	% compute electric field
	Ef=-Jo*rf*cosh(c(:,1).*lambda)+Vo*lambda*sinh(c(:,1).*lambda);
	Ey=zeros(length(Ef),1);
	Ef=[Ef,Ey];
endfunction

% Wrapper function to call sourcefield
function [O,E,In,Ip]=SourceField(c, m, cp, cn, R, Rg, ll, ur,CMD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
%	m: 	Mesh for computing the	[cm]:[cm]
%		total current
%	cp:	coordinates and radius	[cm],[cm],[cm]
%		p-contact x,y,r
%	cn:	coordinates and radius	[cm],[cm],[cm]
%		n-contact x,y,r
%	R:	sheet resistance	[Ohm]
%	Rg:	Junction resistance	[Ohm cm2]
%	ll:	lower left corner	[cm],[cm]
%	ur:	upper right corner	[cm],[cm]	
%	CMD:	Computation Commands
%		for sourcefield
% OUTPUT:
%	O:	generic output data	[?]
%	E:	Electric Field @m	[V/cm]:[V/cm]
%		in x&y direction
%	In:	current through 	[A]
%		n-contact
%	Ip:	current through 	[A]
%		p-contact
	% write a sourcefield input file
	global SFErr;
	fn="dummy_file";
	f=fopen(fn, "w");
		
	fprintf(f,"x1 %e\ny1 %e\n", ll(1),ll(2));
	fprintf(f,"x2 %e\ny2 %e\n", ur(1),ur(2));
	
	% here you can imrpove the speed alot, or not...
	fprintf(f,"err %e\n", SFErr);
	fprintf(f,"rp %e\n", cp(3));
	fprintf(f,"rn %e\n", cn(3));
	fprintf(f,"rg %e\n", Rg);
	fprintf(f,"Rsq %e\n", R);	
	fprintf(f,"cxp %e\ncyp %e\n", cp(1),cp(2));
	fprintf(f,"cxn %e\ncyn %e\n", cn(1),cn(2));
	fprintf(f,"mesh dummy_mesh\n");
	fprintf(f,"file dummy_out\n");
	for i=1:length(CMD(:,1))
		fprintf(f,"%s\n",CMD(i,:));
	endfor
	
	fclose(f);
	
	  
	% write input mesh to dummy_file
	WriteMesh(c,"dummy_mesh");
	WriteMesh(m,"current_mesh");
	
	% Run Sourcefield, dump stdout to a dile
	system("./SourceField dummy_file > dummy_stdout",0);
	% Load the voltages
	O=load("dummy_out");
	E=load("Ex0.dat");
	
	% From this file filter out the actual currents (not always 1A)
	[s,In]=system("egrep -e \'In [^\\ ]+\' -o dummy_stdout|egrep -e \'[^In\\ ].*\' -o");
	[s,Ip]=system("egrep -e \'Ip [^\\ ]+\' -o dummy_stdout|egrep -e \'[^Ip\\ ].*\' -o");
	In=str2num(In);
	Ip=str2num(Ip);
	% exit(1)
	% clean up
	system("rm dummy_file",0);
	system("rm dummy_stdout",0);
	system("rm dummy_mesh",0);
	system("rm current_mesh",0);
	system("rm dummy_out",0);
	system("rm Ex0.dat",0);
endfunction


% to compute the external current induced by a shunt we must integrate the current along the contact
% To do that efficiently we have to generate an appropriate mesh. Here we make an estimate of the current
% density along the contact based on a simple point source in a 2D plane. As the point source sends out 
% current uniformly in all directions it is best to make the mesh equidistant in angle. Alng the contact 
% we then get a mesh that follows a tangens function:
function y=GenMesh(N, cn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	N: 	Number of points
%	cn: 	coordinates and radius	[cm],[cm],[cm]
%		n-contact x,y,r
% OUTPUT:
%	y:	y coordinates of mesh	[cm]
	global x1;
	global y1;
	global y2; 
	
	a=-pi/2+pi/N:pi/N:pi/2-pi/N;		% equidistant mesh in angle
	y=tan(a).*(cn(1,1)-x1).+cn(1,2);	% Corresponmding y coordinates
	y=y';
	ii=(y<y2);				% select the y coordinates within range
	y=y(ii);
	ii=(y>y1);
	y=y(ii);
	y=[y1;y;y2];				% add that start and end points
endfunction


% Wrapper function to call sourcefieldG(
function [V,I,In,Ip]=SourceFieldPotential(c, cp, cn, R, Rg, ll, ur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
%	cp:	coordinates and radius	[cm],[cm],[cm]
%		p-contact x,y,r
%	cn:	coordinates and radius	[cm],[cm],[cm]
%		n-contact x,y,r
%	R:	sheet resistance	[Ohm]
%	Rg:	Junction resistance	[Ohm cm2]
%	ll:	lower left corner	[cm],[cm]
%	ur:	upper right corner	[cm],[cm]
% OUTPUT:
%	V:	Voltage @ c		[V]
%	I:	Total external current	[A]
%	In:	current through 	[A]
%		n-contact
%	Ip:	current through 	[A]
%		p-contact
	global x1;
	global x2;
	global y1;
	global y2;
	global Nnumint;
	m=GenMesh(Nnumint, cn);
	m=[ones(length(m),1).*x1,m];
	CMD=["bess_fast";"bess_potential";"mesh current_mesh";"file Ex0.dat";"bess_field";];	
	[V,E0, In, Ip]=SourceField(c, m, cp, cn, R, Rg, ll, ur,CMD);
	V=V(:,3);
	% integrate the external current
	I0=E0(:,3)./R;
	dx=m(2:end,2).-m(1:end-1,2);
	ya=(I0(2:end)+I0(1:end-1))/2;
	ya.*=dx;
	I=-sum(ya);
endfunction

function [E,I,In,Ip]=SourceFieldEField(c, cp, cn, R, Rg, ll, ur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
%	cp:	coordinates and radius	[cm],[cm],[cm]
%		p-contact x,y,r
%	cn:	coordinates and radius	[cm],[cm],[cm]
%		n-contact x,y,r
%	R:	sheet resistance	[Ohm]
%	Rg:	Junction resistance	[Ohm cm2]
%	ll:	lower left corner	[cm],[cm]
%	ur:	upper right corner	[cm],[cm]
% OUTPUT:
%	E:	Electric field @ c	[V/cm]:[V/cm]
%		in x&y direction
%	I:	Total external current	[A]
%	In:	current through 	[A]
%		n-contact
%	Ip:	current through 	[A]
%		p-contact	
	global x1;
	global x2;
	global y1;
	global y2;
	global Nnumint;
	m=GenMesh(Nnumint, cn);
	m=[ones(length(m),1).*x1,m];
	CMD=["bess_fast";"bess_field";"mesh current_mesh";"file Ex0.dat";"bess_field";];
	[E,E0, In, Ip]=SourceField(c, m, cp, cn, R, Rg, ll, ur,CMD);
	E=E(:,3:4);
	I0=E0(:,3)./R;
	dx=m(2:end,2).-m(1:end-1,2);
	ya=(I0(2:end)+I0(1:end-1))/2;
	ya.*=dx;
	I=-sum(ya);
endfunction
	

% Computes the potential acording to the sourcefield program
% potentials normalized to 1A current through the shunt, external current is returned in Iext
function [V, Iext]=ShuntPotential(c, sh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
%	sh:	coordinates and radius	[cm],[cm],[cm]
%		of shunt x,y,r
% OUTPUT:
%	V:	Voltage @ In=Ip=1A	[V]
%	Iext:	Total external current	[A]
	global Rf;
	global lin_sys;
	global x1;
	global x2;
	global y1;
	global y2;
	Rj=lin_sys(1,4);
	% compute the potentials due to the front and back contact seperately
	% To do that we must divide the junction resistance to the respecting sheet resistance after their ratio
	%                            anti shunt                 shunt                                      image plane, plane
	[Vf,I, In, Ip]=SourceFieldPotential(c, [x1-sh(1),sh(2),sh(3)], [sh(1),sh(2),sh(3)], Rf, Rj,[x1-x2,y1],[x2,y2]);	
	% Normalize to the current through the shunt
	Vf./=Ip;	
	Iext=I./Ip;
	V=Vf;	
endfunction

function Ef=ShuntEField(c, sh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	c: 	columar coordinates x:y	[cm]:[cm]
%	sh:	coordinates and radius	[cm],[cm],[cm]
%		of shunt x,y,r
% OUTPUT:
%	E:	Electric field @ 	[V/cm]:[V/cm]
%		In=Ip=1A in x&y
	global Rf;
	global Rb;
	global lin_sys;
	global x1;
	global x2;
	global y1;
	global y2;
	Rj=lin_sys(1,4);
	for i=1:length(sh(:,1))	
		% compute the potentials due to the front and back contact seperately
		% To do that we must divide the junction resistance to the respecting sheet resistance after their ratio
		%                            anti shunt                 shunt                                      image plane, plane
		[Ef,I, In, Ip]=SourceFieldEField(c, [x1-sh(i,1),sh(i,2),sh(i,3)], [sh(i,1),sh(i,2),sh(i,3)], Rf, Rj,[x1-x2,y1],[x2,y2]);		
		
	
		% Normalize to the current
		Ef./=Ip;		
	endfor
endfunction
% 

% Our various solutions need to be superimposed to get the final solution
% for that I need the current for each solution
% this routine solves the linear system which solves for:
% base-potential gets the diode current
% voltage at each shunt coordinate divided by the shunt resistance is the current through that shunt

function [Vsh, Ishext]=Vshunts()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	Initialized global variable lin_sys
% OUTPUT:
%	All currents and voltages assume a shunt current of 
% 	1A for all shunts
%	Vsh:	Voltage matrix:		[[V]::[V];;[V]::[V]]
%		Vij=Voltage in shunt j 
%		due to shunt i
%	Ishext:	External current for 	[A]
%		shunt i

	global lin_sys;
	
	% number of superimposed potential functions
	N=length(lin_sys(:,1))-1;
	
	Vsh=zeros(N,N);
	Ishext=zeros(N,1);
	for i=1:N
		cc=lin_sys(2:end,1:2);
		[Vshunts, Iext]=ShuntPotential(cc, lin_sys(i+1,1:3));
		Vsh(i,:)=Vshunts';
		Ishext(i)=Iext;
	endfor	
endfunction

% Solves all currents through each shunt given the total current through the diode, Id
function [I,Iext]=ShuntCurr(Vsh, Ishext, Id)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	Vsh:	Voltage matrix:		[[V]::[V];;[V]::[V]]
%		Vij=Voltage in shunt j 
%		due to shunt i
%	Ishext:	External current for 	[A]
%		shunt i
%	Id:	total diode current 	[A]
% OUTPUT:
%	I:	Currents through each 	[A]
%		solution
%	Iext:	External current for 	[A]
%		each solution		[A]
	global lin_sys;
	
	% number of superimposed potential functions
	N=length(lin_sys(:,1));
	
	% b is the right hand side vector and contains all sources
	b=zeros(N, 1).-lin_sys(1,5);
	b(2:end)=(b(2:end).+lin_sys(2:end,5));
	b(1,1)=Id;
	
	A=zeros(N,N);
	% base-potential gets the diode current
	A(1,:)=zeros(1,N);
	A(1,1)=1;
	for i=2:N
		cc=lin_sys(i,1:2);
		Vbase=BasePotential(cc);
		Vshunts=Vsh(i-1,:);
		A(i,:)=[Vbase,Vshunts];
		% A(i,:)./=(lin_sys(i,4))/(pi*lin_sys(i,3)^2);
		A(i,i)-=(lin_sys(i,4))/(pi*lin_sys(i,3)^2);
	endfor	
	ii=(abs(A)<1e-16);
	A(ii)=0;
	I=A\b;
	Iext=I;
	Iext(2:end).*=Ishext;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the voltage and power for each part solution
function V=SolV(I)
% INPUT:
%	Initialized global variable lin_sys
%	current through each part solution
% OUTPUT:
%	Bias voltage for each part solution
	global lin_sys;
	V=I.*lin_sys(:,4)./(pi*lin_sys(:,3).^2).+lin_sys(:,5);
endfunction


% Sets the junction linearization and iteratively solves for each (non-linear) shunt the currents
function [I,Iext]=NLShuntCurr(Id, e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	Id:	total diode current 	[A]
%	e:	Maximum relative error	[-]
% OUTPUT:
%	I:	Currents through each 	[A]
%		solution
%	Iext:	External current for 	[A]
%		each solution		[A]
	global lin_sys;
	global VJR={};
	
	% set junction working point
	[R,Vc]=LinVI_I(Id/(pi*lin_sys(1,3)^2), VJR{1});
	lin_sys(1,4:5)=[R,Vc];
	[Vsh, Ishext]=Vshunts();
	
	% The shunt working points are unaltered for the moment
	[I,Iext]=ShuntCurr(Vsh, Ishext, Id);
	N=length(I);
	Io=I.-2*e;
	% Iterate to the solution
	while (sqrt(sum((I.-Io).^2)/length(Io))>e)
		Io=I;
		for i=2:N
			[R,Vc]=LinVI_I(I(i)/(pi*lin_sys(i,3)^2), VJR{i});
			lin_sys(i,4:5)=[R,Vc];
		endfor	
		[I,Iext]=ShuntCurr(Vsh, Ishext, Id);
	endwhile
	for i=2:N
		[R,Vc]=LinVI_I(I(i)/(pi*lin_sys(i,3)^2), VJR{i});
		lin_sys(i,4:5)=[R,Vc];
	endfor	
endfunction

% computes the potential by superimposing all potential functions
% V is the junction voltage, Vf the front contact voltage and Vb the back contact voltage
function V=Potential(c,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	c: 	columar coordinates x:y	[cm]:[cm]
%	I:	Currents through each 	[A]
%		solution
% OUTPUT:
%	V:	Voltage @ c		[V]
	global lin_sys;
	% number of superimposed potential functions
	N=length(I);
	
	% add the potential functions, scaled to their respective currents
	V=BasePotential(c);
	V=I(1).*V.+lin_sys(1,5); % add the current independent Vc
	for i=2:N
		dV=ShuntPotential(c, lin_sys(i,:));
		V+=dV.*I(i);
	endfor
	V=[c,V];
endfunction

function Ef=EField(c,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	c: 	columar coordinates x:y	[cm]:[cm]
%	I:	Currents through each 	[A]
%		solution
% OUTPUT:
%	Ef:	Electric field @ c	[V/cm]:[V/cm]
%		in x&y
	global lin_sys;
	% number of superimposed potential functions
	N=length(I);
	
	% add the potential functions, scaled to their respective currents
	[V,Ef]=BasePotential(c);
	Ef=I(1).*Ef;
	for i=2:N
		dEf=ShuntEField(c, lin_sys(i,:));
		Ef+=dEf.*I(i);
	endfor
	Ef=[c,Ef];
endfunction

% Computes the external voltage, i.e. voltage difference between the contacts with the front and back electrodes
function V=ExtVoltage(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	I:	Currents through each 	[A]
%		solution
% OUTPUT:
%	V:	External voltage	[V]
	global x1;
	global x2;
	global y1;
	global y2;
	global t_extv;
	global n_extv;
	t1=time();
	n_extv++;
	% solve the potentials
	V=Potential([x1,(y1+y2)/2],I);
	V=V(3);
	t_extv+=(time()-t1);
	
endfunction

function [V,Ef,Jf, Jj, Pf,Pj]=VEJP(c,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	c: 	columar coordinates x:y	[cm]:[cm]
%	I:	Currents through each 	[A]
%		solution
% OUTPUT:
%	V:	potential @ c		[V]
%	Ef:	Electric field @ c	[V/cm]:[V/cm]
%		in x&y
%	Jf:	Current density @ c	[A/cm]:[A/cm]
%		in electrode x&y
%	Jj:	Current density @ c	[A/cm2]
%		through diode
%	Pf:	Power density @ c	[W/cm2]
%		in electrode
%	Pj:	Power densita @ c	[W/cm2]
%		in diode
	global lin_sys;
	global Rf;
	% number of superimposed potential functions
	N=length(I);
	
	% add the potential functions, scaled to their respective currents
	[V,Ef]=BasePotential(c);
	V=I(1).*V.+lin_sys(1,5); % add the current independent Vc
	Ef=I(1).*Ef;
	for i=2:N
		dEf=ShuntEField(c, lin_sys(i,:));
		Ef+=dEf.*I(i);
		dV=ShuntPotential(c, lin_sys(i,:));
		V+=dV.*I(i);
	endfor
	Jf=-Ef./Rf;
	Jj=(V.-lin_sys(1,5))./lin_sys(1,4);
	
	for i=2:N
		d2=(c(:,1).-lin_sys(i,1)).^2+(c(:,2).-lin_sys(i,2)).^2;
		ii=(d2<lin_sys(i,3)^2);
		Jj(ii)=I(i)/(pi*lin_sys(i,3)^2);
	endfor
	
	Pf=-Jf.*Ef;
	Pf=sqrt(Pf(:,1).^2+Pf(:,2).^2);
	Pj=Jj.*V;
	
	V=[c,V];
	Ef=[c,Ef];	
	Jf=[c,Jf];
	Jj=[c,Jj];
	Pf=[c,Pf];
	Pj=[c,Pj];
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Non-Linear iterative solvers
% Solves the device for a total injected current
function [I, Ve, Iext]=SetCurrent(Ii, e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	Ii:	total injected current 	[A]
%	e:	Maximum relative error	[-]
% OUTPUT:
%	I:	Currents through each 	[A]
%		solution
%	Ve:	External voltage	[V]
%	Iext:	External current for 	[A]
%		each solution
	global Nbisect;
	Id1(1)=Ii;
	[I,Id1]=NLShuntCurr(Id1(1), e);	
	E1=Ii-sum(Id1);	
	E1=sum(Id1)-Ii;
	
	Id2(1)=Ii-sum(Id1(2:end));
	[I,Id2]=NLShuntCurr(Id2(1), e);
	E2=sum(Id2)-Ii;
	
	
		
	printf("Solve for a set current of %e\n",Ii);
	printf("Iter\tIdmax\t\tId\t\tIdmin\t\tEmax\t\tE\t\tEmin\n");
	if (E2>E1)
		% swap if needed
		E=E2;
		Id=Id2;
		E2=E1;
		Id2=Id1;
		E1=E;
		Id1=Id;
	endif
	while (E1<0)
		Id1(1)=Ii+sum(Id1(2:end));
		[I,Id1]=NLShuntCurr(Id1(1), e);
		E1=sum(Id1)-Ii;
	endwhile
	
	while (E2>0)
		Id2(1)=Ii-sum(Id2(2:end));
		[I,Id2]=NLShuntCurr(Id2(1), e);
		E2=sum(Id2)-Ii;
	endwhile
	
	
	MaxR=1;
	aEmax=abs(E2);
	aEmin=abs(E1);
	if(aEmin/aEmax>MaxR)
		aEmin=MaxR*aEmax;
	else if (aEmax/aEmin>MaxR)
		aEmin=MaxR*aEmax;
	endif
	endif
	
	Id=(aEmax.*Id1+aEmin.*Id2)./(aEmax.+aEmin);
	% Now we can bisect
	% Id=(-E2*Id1+E1*Id2)/(E1-E2);
	% Id=(Id1+Id2)./2; % this is the most rubust one when the junction resistance becomes very high
	[I,Id]=NLShuntCurr(Id(1), e);
	E=sum(Id)-Ii;	
	iter=1;	
	while (abs(E)>e)
		printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\n", iter,Id1(1),Id(1),Id2(1), E1, E, E2);
		if (E<0)
			E2=E;
			Id2=Id;
		else
			E1=E;
			Id1=Id;
		endif
		aEmax=abs(E2);
		aEmin=abs(E1);
		if(aEmin/aEmax>MaxR)
			aEmin=MaxR*aEmax;
		else if (aEmax/aEmin>MaxR)
			aEmin=MaxR*aEmax;
		endif
		endif
		if (rem(iter,Nbisect))
			Id=(aEmax.*Id1+aEmin.*Id2)./(aEmax.+aEmin);
		else
			Id=(Id1+Id2)./2;
		endif
		[I,Id]=NLShuntCurr(Id(1), e);
		Eo=E;
		E=sum(Id)-Ii;
		MaxR=1+log(max([abs(Eo/E),1]));
		iter++;
	endwhile
	printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\n", iter,Id1(1),Id(1),Id2(1), E1, E, E2);
	Ve=ExtVoltage(I);
	Iext=Id;
endfunction

% Solves the device for an applied voltage
function [I,Iext]=SetVoltage(V, e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	V:	External voltage	[V]
%	e:	Maximum relative error	[-]
% OUTPUT:
%	I:	Currents through each 	[A]
%		solution
%	Iext:	External current for 	[A]
%		each solution
	global x1;
	global x2;
	global y1;
	global y2;
	global VJR={};
	global Nbisect;
	A=(x2-x1)*(y2-y1);
	
	V2=V;
	[Rj,Vc]=LinVI_V(V2,VJR{1});
	Id=A*(V2-Vc)/Rj;
	[I,Iext]=NLShuntCurr(Id, e);
	Ve=ExtVoltage(I);	
	E2=V-Ve;
	if (abs(E2)<e)
		return;
	endif
	
	V1=V+E2;
	[Rj,Vc]=LinVI_V(V1,VJR{1});
	Id=A*(V1-Vc)/Rj;
	[I,Iext]=NLShuntCurr(Id, e);
	Ve=ExtVoltage(I);	
	E1=V-Ve;
	if (abs(E1)<e)
		return;
	endif
	
	printf("Solve for a set voltage of %e\n", V);
	printf("Iter\tVmax\t\tV\t\tVmin\t\tEmax\t\tE\t\tEmin\n");
	if (E2>E1)
		% swap if needed
		E=E2;
		Vn=V2;
		E2=E1;
		V2=V1;
		E1=E;
		V1=Vn;
	endif
	while (E2>0)
		V2+=(V2-V1);
		[Rj,Vc]=LinVI_V(V2,VJR{1});
		Id=A*(V2-Vc)/Rj;
		[I,Iext]=NLShuntCurr(Id, e);
		Ve=ExtVoltage(I);	
		E2=V-Ve;
		if (abs(E2)<e)
			return;
		endif
	endwhile
		
	% Now we can bisect
	Vn=(-E2*V1+E1*V2)/(E1-E2);
	[Rj,Vc]=LinVI_V(Vn,VJR{1});
	Id=A*(Vn-Vc)/Rj;
	[I,Iext]=NLShuntCurr(Id, e);
	Ve=ExtVoltage(I);		
	E=V-Ve;
	
	iter=1;	
	while(abs(E)>e)
		printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\n", iter,V2,Vn,V1, E2, E, E1);
		if (E<0)
			E2=E;
			V2=Vn;
		else
			E1=E;
			V1=Vn;
		endif
		if (rem(iter,Nbisect))
			Vn=(-E2*V1+E1*V2)/(E1-E2);
		else
			Vn=(V1+V2)/2;
		endif
		[Rj,Vc]=LinVI_V(Vn,VJR{1});
		Id=A*(Vn-Vc)/Rj;
		[I,Iext]=NLShuntCurr(Id, e);
		Ve=ExtVoltage(I);		
		E=V-Ve;
		iter++;
	endwhile
	printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\n", iter,V2,Vn,V1, E2, E, E1);
endfunction


% Finds the solution for the maximum powerpoint
function [Impp,Vmpp, Iextmpp]=FindMPP(e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	e:	Maximum relative error	[-]
% OUTPUT:
%	Impp:	Currents through each 	[A]
%		solution
%	Vmpp:	External voltage	[V]
%	Iextmpp:External current for 	[A]
%		each solution
	global x1;
	global x2;
	global y1;
	global y2;
	global VJR={};
	
	printf("Find maximum powerpoint\n");
	A=(x2-x1)*(y2-y1);
	P=zeros(5,1);
	Id=zeros(5,1);
	
	% find mpp of the diode curve
	% at Isc the power is 0
	[I,Iext]=SetVoltage(0, e);
	Id(1)=I(1);
	
	[I,V,Iext]=SetCurrent(sum(Iext)/2, e);
	Id(3)=I(1);
	P(3)=-sum(Iext)*V
	
	Id(5)=-sum(I(2:end));
	iter=0;
	printf("Iter\tP1\t\tP\t\tP2\n");
	while (Id(5)-Id(1)>2*e)
		printf("%i\t%e\t%e\t%e\n", iter,P(1), P(3), P(5));
		Id(2)=(Id(3)+Id(1))/2;
		[I,Iext]=NLShuntCurr(Id(2), e);
		V=ExtVoltage(I);
		P(2)=-sum(Iext)*V;
		
		Id(4)=(Id(3)+Id(5))/2;
		[I,Iext]=NLShuntCurr(Id(4), e);
		V=ExtVoltage(I);
		P(4)=-sum(Iext)*V;
		
		[Pm,im]=max(P);
		P([1,3,5])=P([im-1,im,im+1]);
		Id([1,3,5])=Id([im-1,im,im+1]);	
		iter++;	
	endwhile
	[Impp, Iextmpp]=NLShuntCurr(Id(3), e);
	Vmpp=ExtVoltage(Impp);	
endfunction


function SetDiode(VJ)
	global x1;
	global x2;
	global y1;
	global y2;
	global VJR;
	global lin_sys;
	
	% From the JV data we determine the resistance versus voltage/current density
	R=Diff(VJ(:,1))./Diff(VJ(:,2));
	% add this JV to the VJR list:
	VJR{1}=[VJ,R];
	% Linearize the diode
	[R,Vc]=LinVI_V(0,VJR{1});
	% add to our linear system matrix
	lin_sys(1,:)=[0,0,sqrt((x2-x1)*(y2-y1)/pi),R,Vc]; % diode at 0,0 (this is needless info, just for indexing it is convenient
endfunction

function AddPointSource(VJ,x,y,r)
	global x1;
	global x2;
	global y1;
	global y2;
	global VJR;
	global lin_sys;
	
	if (length(lin_sys(:,1))==0)
		printf("Error, no diode defined yet. Please use the \"SetDiode\" command first\n")
		exit;
	endif
	
	R=Diff(VJ(:,1))./Diff(VJ(:,2));
	VJR{end+1}=[VJ,R];
	[R,Vc]=LinVI_V(0,VJR{end});
	lin_sys(end+1,:)=[x,y, r,R,Vc];
	printf("Created Point source with index %d\n",length(lin_sys(:,1))-1);
endfunction


function SetPointSource(index,VJ,x,y,r)
	global x1;
	global x2;
	global y1;
	global y2;
	global VJR;
	global lin_sys;
	
	if (length(lin_sys(:,1))<index+1)
		printf("Error, this JV characteristic does not exist yet\n")
		exit;
	endif
	
	R=Diff(VJ(:,1))./Diff(VJ(:,2));
	VJR{index+1}=[VJ,R];
	[R,Vc]=LinVI_V(0,VJR{index+1});
	lin_sys(index+1,:)=[x,y, r,R,Vc];
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Useless Utillity Functions
% A progress bar, you never know when you need one.
function ProgressScale(len, tics)
	tic=floor(len/tics);
	for i=0:len
		if(rem (i, tic) == 0)
			printf("|%-3.0f", 100*i/len);
		else if (rem(i,tic)>3)
			printf(" ");
		endif
		endif
	endfor

	printf("\n");
endfunction
function pc=ProgressBar(pcn, pco, len, tics)	
	% input pcn: new percentage complete 
	% input pco: old percentage complete 
	% input len: length of the progress bar
	% output (new) old percentage	
	
	% This function is supposed to be used in a loop of known length
	% Usage is best explained in an example:
	% pco=-1;
	% for i=1:N
	%     pcn=100*i/N;
	%     pco=ProgressBar(pcn, pco, len, tics) %len is the length of the progress bar and tics the number of tics
	%     ...
	% end
	%
	
	% round to integer and norm to length
	if (pco<0)
		ProgressScale(len, tics);
	endif
	
	pc_n=floor(pcn*len/100);	
	pc=pco;		
	pco=floor(pco*len/100);
	if (pco==len)
		return;
	endif
	tic=floor(len/tics);
	if(pc_n>pco)
		if(rem (pc_n, tic) == 0)
			if (pc_n==0)
				printf("|");
			else
				printf("\b\b\b\b");
				printf("|");
			endif
		else
			printf("\b\b\b\b");
			printf("=");
		endif
		pc=pcn;
		printf("%3.0f%%",pc);
	else
		if (floor(pcn)>floor(pc))
			printf("\b\b\b\b");
			printf("%3.0f%%",pcn);		
		endif
	endif
	pc=pcn;
	if (pc_n==len)
		printf("\n");
	endif
endfunction
% save data, you can use the built-in save but why bother?
function SaveXY(file, data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%	file: 	filename
%	data: 	Matrix with data to be written

	fid=fopen(file,"w");	
	for i=1:length(data(:,1))
		for j=1:length(data(1,:))
			fprintf(fid,"%.12e\t",data(i,j));			
		endfor	
		fprintf(fid, "\n");
	endfor	
	fclose(fid);
endfunction 
