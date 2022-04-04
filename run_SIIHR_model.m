%Script for running and plotting deterministic and stochastic SIIHR ODE models
clear all
clf
%Define model parameters as a structure
para = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
%Define initial conditions as a structure
ICs = struct('S',para.N-1,'Ia',0,'Is',1,'H',0,'D',0,'R',0);

%Define time to run model for
mintime=0;
maxtime = 500;

%Run model by calling function ODE_SIIHRmodel.m
[Classes] = ODE_SIIHR_model(para,ICs,mintime,maxtime);

%Q1 a)
%Plot dynamics

Cmap=[0,128,255;255,128,0;204,0,0;255,102,178;102,0,204;0 0 0;0,153,0]./255;
colororder(Cmap)
figure(1)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',3)
plot(Classes.t,Classes.S,Classes.t,Classes.Ia,Classes.t,Classes.Is,Classes.t,Classes.Ia+Classes.Is,Classes.t,Classes.H,Classes.t,Classes.D,Classes.t,Classes.R)
xlabel('Time (days)')
ylabel('Number of people')
legend('S','I_{asymptomatic}','I_{symptomatic}','I_{total}','H','D','R')
title('SI_aI_sIHDR Dynamics with given parameters')
ylim([0 70e6])
xlim([0 350])

% Compute epidemic quantities
% final size
Final_size=round(Classes.R(end)+ Classes.D(end))%use round() to get to the nearest whole number

Death_size=round(Classes.D(end))
%Compute expected outbreak duration
Last_Inft = find(Classes.Ia+Classes.Is>1,1,'last'); %finds the last remaining infection
Duration = ceil(Classes.t(Last_Inft))

%peak of outbreak
 peak=find(Classes.Ia+Classes.Is==max(Classes.Ia+Classes.Is))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Running tau leap model for SIIHDR model

%Define the time step in tau leaping algorithm
timestep = 1;

%Run stochastic model using tau leap algorithm with same ICS and params as
%deterministic model
tic  %timing how long it takes matlab to do one iteration using tau leap
[Classes_TL] = Tauleap_SIIHDR_model(para,ICs,maxtime,timestep);
toc

%Plot one realisation of infections from tau leap algorithm with ODE
%infections on same graph
figure(2)
plot(Classes_TL.t, (Classes_TL.Ia +Classes_TL.Is),'r',Classes.t,Classes.Ia+Classes.Is,'--k')
xlabel('Time (days)')
ylabel('Infections (I)')
legend('Tau Leap','ODE')
title('A sample path of the stochastic model')



%%
%Plot 500 realisations of tau leap model

%Plot your existing stochastic tau leap realisation on the graph now
figure(3)
plot(Classes_TL.t,Classes_TL.Ia,'r')
hold on

%Use a for loop to run each tau leap simulation, store metrics then plot the
%dynamics
NoRuns=500;

for i=2:NoRuns
    
    [Classes_TL] = Tauleap_SIIHDR_model(para,ICs,maxtime,timestep);
    
    final_sizeT(i)=Classes_TL.R(end)+ Classes_TL.D(end);
    %plot next 499 realisations
    h3=plot(Classes_TL.t,Classes_TL.Ia,'color','r');
   
end


%Add the ODE solution back on top of the stochastic solutions
h4=plot(Classes.t,Classes.Ia,'--k','linewidth',3);
axis([0 maxtime 0 Inf])
xlabel('Time (days)')
ylabel('Infections (I)')
xlim([0 maxtime])
ylim([0 inf])
title("500 Realisations using Tau Leap Model")
legend([h3 h4],'Tau leap','ODE')


%%

%%Plot histograms of final size
figure(3)
clf
histogram(final_sizeT,'numbins',50,'normalization','probability')
xlabel('Final size')
ylabel('Probability')
hold on
ymax=get(gca,'ylim');
plot([Final_size Final_size], [0 ymax(2)])
legend('500 tau leap realisations','ODE')
title('Distribution of final size')



