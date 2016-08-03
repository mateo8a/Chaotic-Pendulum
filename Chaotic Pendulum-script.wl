(* ::Package:: *)

Final Project


Scientific Computation


Mateo Ochoa


Functions


The Runge-Kutta 4 method is given by:


RKODE4[G_, x0_,f0_, dx_, Ns_] := Module[{k1,k2,k3,k4,retvals,index,xn,fn},
retvals = Table[0.0,{i,1,Ns}];
xn = x0;
fn = f0;
retvals[[1]] = {x0,f0};
For[index = 2, index <= Ns, index=index+1,
k1 = dx G[xn,fn];
k2 = dx G[xn + dx/2, fn + k1/2];
k3=dx G[xn+dx/2,fn+k2/2];
k4=dx G[xn+dx,fn+k3];
fn = fn + 1/3(k1/2 +k2+k3+k4/2);
xn = xn + dx;
retvals[[index]] = {xn,fn};
];
Return[retvals];
]


We limit the range of the angular displacement values to live between -Pi and Pi:


LimitAngleRange[data_]:=Module[{moddata,i},
moddata=data;
For[i=1,i<=Length[data],i=i+1,
While[moddata[[i,2,1]]> Pi,
moddata[[i,2,1]]=moddata[[i,2,1]]-2Pi];
];
For[i=1,i<=Length[data],i=i+1,
While[moddata[[i,2,1]]<-Pi,
moddata[[i,2,1]]=moddata[[i,2,1]]+2Pi];
];
Return[moddata];
]


We can plot the Phase Space using:


PhaseSpacePlot[moddata_]:=Module[{i,j,k,plot,anglevsw},
anglevsw = Table[{moddata[[j,2,1]],moddata[[j,2,2]]},{j,1,Length[moddata]}];
plot=ListPlot[anglevsw,PlotStyle->PointSize[Tiny],AxesLabel->{"\[Theta]","\[Omega]"}];
Return[plot];
]


For easier visualization of the dynamics, we recur to plot the Poincare Section:


PoincareSection[moddata_]:=Module[{eps,i,k,kk,n,list,list2,plot,min,max},
list=Table[0.0,{i,1,Length[moddata]}];
eps=0.05;
n=(Round[moddata[[1,2,3]],2Pi]+2Pi)/(2Pi);
k=1;
min=0;
max=0;

For[i=1,i<=Length[moddata],i=i+1,
If[moddata[[i,2,3]]<2n Pi+eps && moddata[[i,2,3]]>2n Pi -eps,
list[[k]]={moddata[[i,2,1]],moddata[[i,2,2]]};
n=n+1;
k=k+1;
];
];

For[i=1,i<=Length[moddata],i=i+1,
If[moddata[[i,2,2]]<min,min=moddata[[i,2,2]]];
If[moddata[[i,2,2]]>max,max=moddata[[i,2,2]]];
];

kk=k-1;
list2=Table[0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
list2[[i]]=list[[i]]
];
plot=ListPlot[list2,PlotStyle->PointSize[0.01],PlotRange->{All,{min,max}},AxesLabel->{"\[Theta]","\[Omega]"}];
Return[plot]]


In case we only want the data needed to plot a Poincare Section, we can use:


PoincareSectionData[moddata_]:=Module[{eps,i,k,kk,n,list,list2,plot,min,max},
list=Table[0.0,{i,1,Length[moddata]}];
eps=0.05;
n=(Round[moddata[[1,2,3]],2 Pi]+2Pi)/(2 Pi);
k=1;
min=0;
max=0;
For[i=1,i<=Length[moddata],i=i+1,
If[moddata[[i,2,3]]<n 2Pi+eps && moddata[[i,2,3]]> n 2 Pi -eps,
list[[k]]={moddata[[i,2,1]],moddata[[i,2,2]]};
n=n+1;
k=k+1;];
];
kk=k-1;
list2=Table[0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
list2[[i]]=list[[i]]];
Return[list2]]


To get rid of the transient motion and only get the stady-state motion, we pick only the last N steps with the following function:


LastNVals[N_,list_]:=Module[{i,k,list2},
k=1;
list2=Table[0.0,{i,1,N}];
For[i=Length[list]-N+1,i<=Length[list],i=i+1,
list2[[k]]=list[[i]];
k=k+1];
Return[list2];
]


The Fast Fourier Transform:


FFT[invec_] := Module[{k,elist,olist,netlist,eft,oft,n,omega},
If[Length[invec] ==  1, 
Return[invec];
];
n = Length[invec];
omega = 1;
elist = Table[invec[[j]],{j,1,n,2}];
olist = Table[invec[[j]],{j,2,n,2}];
eft = FFT[elist];
oft = FFT[olist];
netlist = Table[0,{j,0,n-1}];
For[k = 0, k <= n/2 - 1, k = k+1,
netlist[[k+1]] = eft[[k+1]] + omega oft[[k+1]];
netlist[[k + n/2+1]] = eft[[k+1]] - omega oft[[k+1]];
omega = omega Exp[2 Pi I/n];
];
Return[netlist];
]


Find the frequency for which the power spectrum peaked (only works with regular motion and if the system has only one driving frequency):


FindSmallest[pspec_]:=Module[{i,max},
For[i=1,i<=Length[pspec],i=i+1,
If[pspec[[i,2]]>1*10^7,max=pspec[[i,1]]];
];
Return[max];
]


To obtain the relevant Power Spectrum data, we filter it and only pick the peaks that are bigger in magnitude than \!\(TraditionalForm\`
\*SuperscriptBox[\(10\), \(8\)]\):


Fk[data_]:=Module[{i,list,listf,k,kk},
list=Table[0,{i,1,Length[data]}];
k=1;
For[i=1,i<=Length[data],i=i+1,
If[(Abs[data[[i]]])^2>10^8,list[[k]]={i-1,data[[i]]};k=k+1]
];

kk=k-1;
listf=Table[0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
listf[[i]]=list[[i]]];
Return[listf]]


Gets the relevant frequencies from the Power Spectrum:


freqHz[data_]:=Module[{i,list,list2,listf,k,kk},
list=Table[0.0,{i,1,Length[data]}];
k=1;
For[i=1,i<=Length[data],i=i+1,
list[[k]]={(-n df/2)+data[[i,1]] df,data[[i,2]]};k=k+1];

k=1;
list2=Table[0.0,{i,1,Length[list]}];
For[i=1,i<=Length[list],i=i+1,
If[list[[i,1]]>=0,list2[[k]]=list[[i]];k=k+1]
];

kk=k-1;
listf=Table[0.0,{i,1,kk}];
For[i=1,i<=kk,i=i+1,
listf[[i]]=list2[[i]]];
Return[listf]]


Gets the data necessary to plot bifurcation diagrams:


BifurcationData[G_,f_]:=Module[{n,dt,Msteps,Nsteps,t,k,listmod,trajectory,psection,i,list,list2},
dt=wd/64;
Nsteps=2^15;
list=Table[0.0,{i,1,0.5Nsteps}];
k=0;
trajectory=RKODE4[G,0.,f,dt,Nsteps];
trajectory=LastNVals[0.5Nsteps,trajectory];
trajectory=LimitAngleRange[trajectory];
psection=PoincareSectionData[trajectory];
For[i=1,i<=Length[psection],i=i+1,
list[[i]]=psection[[i,2]];k=k+1;
];
listmod=Table[0.0,{i,1,k}];
For[i=1,i<=k,i=i+1,
listmod[[i]]=list[[i]]];
Return[listmod];
]


Bifurcation data with varying parameter A:


BifurcationDataA[f_,Ao_,Af_,dA_,B_,wd_]:=Module[{G,n,A,Msteps,Nsteps,t,k,listmod,trajectory,psection,i,list2},
Msteps=Round[(Af-Ao)/dA]+1;
Nsteps=2^15;
list2=Table[0.0,{i,1,Msteps}];
For[n=0,n<=Msteps-1,n=n+1,
A=Ao+n dA;
G[tt_,ff_]:={ff[[2]],A Cos[ff[[3]]]-(1/B) ff[[2]]-Sin[ff[[1]]],wd};
listmod=BifurcationData[G,f];
list2[[n+1]]={A,listmod};
];
Return[list2];
]


Bifurcation data with varying parameter B:


BifurcationDataB[f_,Bo_,Bf_,dB_,A_,wd_]:=Module[{G,B,n,Msteps,Nsteps,t,k,listmod,trajectory,psection,i,list2},
Msteps=Round[(Bf-Bo)/dB]+1;
Nsteps=2^15;
list2=Table[0.0,{i,1,Msteps}];
For[n=0,n<=Msteps-1,n=n+1,
B=Bo+n dB;
G[tt_,ff_]:={ff[[2]],A Cos[ff[[3]]]-(1/B) ff[[2]]-Sin[ff[[1]]],wd};
listmod=BifurcationData[G,f];
list2[[n+1]]={B,listmod};
];
Return[list2];
]


Bifurcation data with varying parameter \!\(TraditionalForm\`
\*SubscriptBox[\(\[Omega]\), \(d\)]\):


BifurcationDatawd[f_,wdo_,wdf_,dwd_,A_,B_]:=Module[{G,wd,n,Msteps,Nsteps,t,k,listmod,trajectory,psection,i,list2},
Msteps=Round[(wdf-wdo)/dwd]+1;
Nsteps=2^16;
list2=Table[0.0,{i,1,Msteps}];
For[n=0,n<=Msteps-1,n=n+1,
wd=wdo+n dwd;
G[tt_,ff_]:={ff[[2]],A Cos[ff[[3]]]-(1/B) ff[[2]]-Sin[ff[[1]]],wd};
listmod=BifurcationData[G,f];
list2[[n+1]]={wd,listmod};
];
Return[list2];
]


Function to plot the bifurcation data:


PlotFurcate[avals_] := Module[{ret,pts,index,entry,jndex},
pts = {PointSize[Tiny]};
For[index = 1, index <= Length[avals], index=index+1,
entry = avals[[index]];
pts = Join[pts,Table[Point[{entry[[1]], entry[[2,j]]}],{j,1,Length[entry[[2]]]}]];
];
ret = Show[Graphics[{PointSize[.01/2],pts}],Axes->True,AxesLabel->{"","\[Omega]"},AxesOrigin->Automatic];
Return[ret];
]


Analysis


The equations for the damped driven nonlinear pendulum is given by:


G[t_,f_] :={f[[2]],A Cos[f[[3]]]-(1/B) f[[2]]-Sin[f[[1]]],wd}


The following parameters can be changed. Their value will determine if the motion is regular or chaotic.


A=0.5;
B=2;
wd=2/3;
f={Pi/2,0,0};
dt=wd/64;
Nsteps=2^17;


We simulate the trajectory with RK4:


trajectory=RKODE4[G,0.,f,dt,Nsteps];
trajectory=LastNVals[0.5Nsteps,trajectory];
trajectory=LimitAngleRange[trajectory];


We plot the trajectory as a function of time:


anglevst = Table[{trajectory[[j,1]],trajectory[[j,2,1]]},{j,1,Length[trajectory]}];
ListPlot[anglevst,AxesLabel->{"t","\[Theta]"},
PlotRange->{{trajectory[[1,1]],trajectory[[Round[0.2Length[trajectory]],1]]},All}]


The phase space of the trajectory:


PhaseSpacePlot[trajectory]


The Poncare Section of the trajectory:


PoincareSection[trajectory]


We use the FFT to get the Power Spectrum of the system:


n=Length[trajectory];
dt=Pi wd/64;
df=1/(n dt);
data=Table[trajectory[[j,2,1]],{j,1,Length[trajectory]}];


FFTtransform=FFT[data];
FFTtransform=RotateRight[FFTtransform,n/2];
pspec=Table[{k df,(Abs[FFTtransform[[(k+n/2)+1]]])^2},{k,-n/2,n/2 -1}];
max=FindSmallest[pspec];


ListLinePlot[pspec,PlotRange->{{0,0.06},All},AxesLabel->{"\!\(\*SubscriptBox[\(f\), \(k\)]\)","|P(\!\(\*SubscriptBox[\(f\), \(k\)]\))\!\(\*SuperscriptBox[\(|\), \(2\)]\)"}]


Show[%190,PlotLabel->"Periodic motion"]


 


We get the relevant frequencies of the power spectrum:


fk=Fk[FFTtransform];
fnc=N[Chop[freqHz[fk]]];
MatrixForm[fnc]


We generate a bifurcation plot by varying A (caution: the computation time for these plots is around 10 hours): 


BifurcationA=BifurcationDataA[f,1,1.5,0.0001,B,wd];
PlotFurcate[BifurcationA]


We generate a bifurcation plot by varying B (caution: the computation time for these plots is around 10 hours): 


BifurcationB=BifurcationDataB[f,1,7,0.0012,A,wd];
PlotFurcate[BifurcationB]


We generate a bifurcation plot by varying \!\(TraditionalForm\`
\*SubscriptBox[\(\[Omega]\), \(d\)]\) (caution: the computation time for these plots is around 10 hours): 


Bifurcationwd=BifurcationDatawd[f,0,1.2,0.00024,A,B];
PlotFurcate[Bifurcationwd]


Show[%23,ImageSize->{560,300},AspectRatio->Full,AxesOrigin->{0,0},AxesLabel->{"\!\(\*SubscriptBox[\(\[Omega]\), \(D\)]\)","\[Omega]"},PlotLabel->"Bifurcation diagram for varying \!\(\*SubscriptBox[\(\[Omega]\), \(D\)]\), A=1.5 and B=2."]


The following parameters were hand-picked because with them, the system will become chaotic: 


A=1.5;
B=2;
wd=2/3;
f={Pi/2,0,0};
dt=wd/64;
Nsteps=2^23;


We simulate the trajectory using RK4:


trajectory=RKODE4[G,0.,f,dt,Nsteps];
trajectory=LastNVals[0.5Nsteps,trajectory];
trajectory=LimitAngleRange[trajectory];


anglevst = Table[{trajectory[[j,1]],trajectory[[j,2,1]]},{j,1,Length[trajectory]}];
ListPlot[anglevst,AxesLabel->{"t","\[Theta]"},
PlotRange->{{trajectory[[1,1]],trajectory[[Round[0.2Length[trajectory]],1]]},All}]


We plot the trajectory as a function of time:


The Poincare Section of the trajectory:


PoincareSection[trajectory]


Show[%28,PlotLabel->"Chaotic motion"]


Some more bifurcation plots:


BifurcationA=BifurcationDataA[f,1,1.5,0.0001,B,wd];
PlotFurcate[BifurcationA]


Show[%87,ImageSize->{460,300},AspectRatio->Full,AxesLabel->{"A","\[Omega]"},PlotLabel->"Bifurcation diagram for varying A, and B=2 and \!\(\*SubscriptBox[\(\[Omega]\), \(D\)]\)=2/3."]
