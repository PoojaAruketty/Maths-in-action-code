%Script for running and plotting SIIHDR ODE with vaccine models
clear all
clf

%Define model parameters as a structure
para = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0,'alpha',0,'ro',0); 

%Define initial conditions as a structure
ICs = struct('S',para.N-1,'Sv',0,'Ia',0,'Is',1,'H',0,'D',0,'R',0);

%Define time to run model for
mintime = 0;
period1=80;
maxtime = period1;

%Run model by calling function ODE_SIIHRmodel.m
[Classes] = ODE_SIIHRv_model(para,ICs,mintime,maxtime);



% Now define new parameter and ICs for next phase
param2=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0.005,'alpha',0.61,'ro',0.005); 
ICs2= struct('S',Classes.S(end),'Sv',Classes.Sv(end),'Ia',Classes.Ia(end),'Is',Classes.Is(end),'H',Classes.H(end),'D',Classes.D(end),'R',Classes.R(end));

%Define new start and end times for period 2 (with restrictions)
mintime = period1;
maxtime = 350;

%Continue simulation following restriction
[Classes2] = ODE_SIIHRv_model(param2,ICs2,mintime,maxtime);

[Classes3] = ODE_SIIHR_model(para,ICs,0,maxtime);

%New period
mintime = 0;
period1=100;
maxtime = period1;

%Run model by calling function ODE_SIIHRmodel.m
[Classest2] = ODE_SIIHRv_model(para,ICs,mintime,maxtime);
% Now define new parameter and ICs for next phase
param2=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0.005,'alpha',0.61,'ro',0.005); 
ICs2= struct('S',Classest2.S(end),'Sv',Classest2.Sv(end),'Ia',Classest2.Ia(end),'Is',Classest2.Is(end),'H',Classest2.H(end),'D',Classest2.D(end),'R',Classest2.R(end));

%Define new start and end times for period 2 (with restrictions)
mintime = period1;
maxtime = 300;

%Continue simulation following restriction
[Classes2t2] = ODE_SIIHRv_model(param2,ICs2,mintime,maxtime);

%3rd period
%New period
mintime = 0;
period1=120;
maxtime = period1;

%Run model by calling function ODE_SIIHRv_model.m
[Classest3] = ODE_SIIHRv_model(para,ICs,mintime,maxtime);
% Now define new parameter and ICs for next phase
param2=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0.005,'alpha',0.61,'ro',0.005); 
ICs2= struct('S',Classest3.S(end),'Sv',Classest3.Sv(end),'Ia',Classest3.Ia(end),'Is',Classest3.Is(end),'H',Classest3.H(end),'D',Classest3.D(end),'R',Classest3.R(end));

%Define new start and end times for period 2 (with restrictions)
mintime = period1;
maxtime = 300;

%Continue simulation following restriction
[Classes2t3] = ODE_SIIHRv_model(param2,ICs2,mintime,maxtime);

%4th period
%New period
mintime = 0;
period1=160;
maxtime = period1;

%Run model by calling function ODE_SIIHRv_model.m
[Classest4] = ODE_SIIHRv_model(para,ICs,mintime,maxtime);
% Now define new parameter and ICs for next phase
param2=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0.005,'alpha',0.61,'ro',0.005); 
ICs2= struct('S',Classest4.S(end),'Sv',Classest4.Sv(end),'Ia',Classest4.Ia(end),'Is',Classest4.Is(end),'H',Classest4.H(end),'D',Classest4.D(end),'R',Classest4.R(end));

%Define new start and end times for period 2 (with restrictions)
mintime = period1;
maxtime = 300;

%Continue simulation following restriction
[Classes2t4] = ODE_SIIHRv_model(param2,ICs2,mintime,maxtime);


%Plot the dynamics for each of the 4 cases
figure(2)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',3)

subplot(4,1,1)
plot(Classes.t,Classes.Ia+Classes.Is,Classes2.t,Classes2.Ia+Classes2.Is)
hold on
plot(Classes3.t,Classes3.Ia+Classes3.Is,'--','Color',[255/255,102/255,178/255])
legend('pre-vaccine','vaccine at 80','no vaccine','Location','northwest')
xlabel('Time (days)')
ylabel('Number of people')
xlim([0 300])
title('Infection dynamics with vaccine introduced at varying times')

subplot(4,1,2)
plot(Classest2.t,Classest2.Ia+Classest2.Is,Classes2t2.t,Classes2t2.Ia+Classes2t2.Is)
hold on
plot(Classes3.t,Classes3.Ia+Classes3.Is,'--','Color',[255/255,102/255,178/255])
legend('pre-vaccine','vaccine at 100','no vaccine','Location','northwest')
xlabel('Time (days)')
ylabel('Number of people')
xlim([0 300])

subplot(4,1,3)
plot(Classest3.t,Classest3.Ia+Classest3.Is,Classes2t3.t,Classes2t3.Ia+Classes2t3.Is)
hold on
plot(Classes3.t,Classes3.Ia+Classes3.Is,'--','Color',[255/255,102/255,178/255])
legend('pre-vaccine','vaccine at 120','no vaccine','Location','northwest')
xlabel('Time (days)')
ylabel('Number of people')
xlim([0 300])

subplot(4,1,4)
plot(Classest4.t,Classest4.Ia+Classest4.Is,Classes2t4.t,Classes2t4.Ia+Classes2t4.Is)
hold on
plot(Classes3.t,Classes3.Ia+Classes3.Is,'--','Color',[255/255,102/255,178/255])
legend('pre-vaccine','vaccine at 160','no vaccine','Location','northwest')
xlabel('Time (days)')
ylabel('Number of people')
xlim([0 300])