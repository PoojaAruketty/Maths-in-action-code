%Script for running and plotting SIIHR ODE models with quarantine

%Define model parameters as a structure
para = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'d_1',0.05,'tau',0.1); 

%Define initial conditions as a structure
ICs = struct('S',para.N-1,'Ia',0,'Is',1,'H',0,'D',0,'R',0,'Q',0);

%Define time to run model for
mintime = 0;
maxtime = 1000;

%Run model by calling function ODE_SIIHRq_model.m
[Classes] = ODE_SIIHRq_model(para,ICs,mintime,maxtime);

%Plot dynamics
Cmap=[0,0,255;255,128,0;204,0,0;255,102,178;102,0,204;0 0 0;0,153,0]./255;
colororder(Cmap)
figure(1)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinelinewidth',3)
plot(Classes.t,Classes.Ia,Classes.t,Classes.Is,Classes.t,Classes.Ia+Classes.Is,Classes.t,Classes.H,Classes.t,Classes.D,Classes.t,Classes.R,Classes.t,Classes.Q)
xlabel('Time (days)')
ylabel('Number of people')
legend('I_{asymptomatic}','I_{symptomatic}','I_{total}','H','D','R','Q')
title('I_aI_sIHDR Dynamics under quarantine')
ylim([0 inf])
xlim([0 500])


% Compute epidemic quantities
% final size
Final_size=round(Classes.R(end)+ Classes.D(end));%use round() to get to the nearest whole number
%Deaths
Death_size=round(Classes.D(end));
%Compute expected outbreak duration
Last_Inft = find(Classes.Ia+Classes.Is>1,1,'last'); %finds the last remaining infection
Duration = ceil(Classes.t(Last_Inft));
%peak
 peak=find(Classes.Ia+Classes.Is==max(Classes.Ia+Classes.Is));


 %Vary d_1 from 0.001 to 0.1
 for t=1:100

        mintime=0;
        maxtime=1200;

      %Define model parameters as a structure
        para2 = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'d_1',t/1000,'tau',1/10); 

      %Define initial conditions as a structure
        ICs = struct('S',para.N-1,'Ia',0,'Is',1,'H',0,'D',0,'R',0,'Q',0);

        %Run the model under quarantine
    [ClassesQ] = ODE_SIIHRq_model(para2,ICs,mintime,maxtime);
    
    %Store final size and duration in arrays
        Final_sizeq1(t)=round(ClassesQ.R(end)+ClassesQ.D(end));
        Last_Inft = find(ClassesQ.Ia+ClassesQ.Is>1,1,'last'); %finds the last remaining infection
        Durationq1(t) = ceil(ClassesQ.t(Last_Inft));
 end

 %
 %Define d1
d1=[1:100]./1000;

%Plot graphs for final size and duration against d1
figure(2)
subplot(1,2,1)
plot(d1,Final_sizeq1)
xlabel('Rate of detection and quarantining (d_1)')
ylabel('Final size')
title('Final size of outbreak while varying d_1')

subplot(1,2,2)
plot(d1,Durationq1)
xlabel('Rate of detection and quarantining (d_1)')
ylabel('Duration')
title('Duration of outbreak while varying d_1')