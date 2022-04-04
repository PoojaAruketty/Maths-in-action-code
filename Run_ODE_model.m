%Script for running and plotting ODE SIR and SIS models

%Script for running and plotting ODE SEIR model for Q1

%Define model parameters as a structure
para = struct('beta',0.3,'gamma',1/10,'N',68000000); 

%Define initial conditions as a structure
ICs = struct('S',para.N-1,'R',0,'I',1);

%Define time to run model for
maxtime = 350;

%Run model by calling function ODE_SIRmodel.m
[Classes] = ODE_SIR_model(para,ICs,maxtime);

%Q1 a)
%Plot dynamics
figure(1)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2.5)
plot(Classes.t,Classes.S,Classes.t,Classes.I,Classes.t,Classes.R)
xlabel('Time (days)')
ylabel('Number of people')
legend('S','I','R')
title('SIR Dynamics with given parameters')
ylim([0 70e6])
xlim([0 inf])


%Compute expected outbreak duration
Last_Inft = find(Classes.I>1,1,'last'); %finds the last remaining infection
Duration = ceil(Classes.t(Last_Inft))% to find the rounded no. of days
Final_size=round(Classes.R(end))%use round() to get to the nearest whole number

 peak=find(Classes.I==max(Classes.I))
 %%
%Define model parameters as a structure
para = struct('beta',0.3,'gamma',1/10,'N',68000000); 

%Define initial conditions as a structure
ICs = struct('S',para.N-1,'R',0,'I',1);

%Define time to run model for
maxtime = 350;

%Run model by calling function ODE_SIRmodel.m
[ClassesS] = ODE_SIS_model(para,ICs,maxtime);

%Q1 a)
%Plot dynamics
figure(2)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',2.5)
plot(ClassesS.t,ClassesS.S,ClassesS.t,ClassesS.I)
xlabel('Time (days)')
ylabel('Number of people')
legend('S','I')
title('SIS Dynamics with given parameters')
ylim([0 70e6])
xlim([0 inf])
