%Modelling lockdown
%Script for running and plotting SIIHR ODE models with lockdown
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
%Define start and end times
period1=90;
mintime=0;
maxtime=period1;

%Rerun with no restriction for 1st period (same param), introduce a "mintime"
[Classes2] = ODE_SIIHR_model(para,ICs,mintime,maxtime);

%Now define new parameter and ICs for next phase
param2=struct('beta1',0.3*0.5,'beta2',0.2*0.5,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
ICs2= struct('S',Classes2.S(end),'Ia',Classes2.Ia(end),'Is',Classes2.Is(end),'H',Classes2.H(end),'D',Classes2.D(end),'R',Classes2.R(end));

%Define new start and end times for period 2 (with restrictions)
mintime = period1;
maxtime = period1+90;

%Continue simulation following restriction
[Classes3] = ODE_SIIHR_model(param2,ICs2,mintime,maxtime);

%Define new initial conditions and parameter for (post restriction) period 3
param3=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
ICs3=struct('S',Classes3.S(end),'Ia',Classes3.Ia(end),'Is',Classes3.Is(end),'H',Classes3.H(end),'D',Classes3.D(end),'R',Classes3.R(end));

%Define start and end times for period 3
mintime = period1+90;
maxtime = 1500;

%Continue simulation following restriction ending
[Classes4] = ODE_SIIHR_model(param3,ICs3,mintime,maxtime);



%Plot the dynamics
figure(2)
set(0,'defaultaxesfontsize',16)
plot(Classes2.t,Classes2.Ia+Classes2.Is,Classes3.t,Classes3.Ia+Classes3.Is,Classes4.t,Classes4.Ia+Classes4.Is)
hold on
plot(Classes.t,Classes.Ia+Classes.Is,'--')
legend('Pre-lockdown','Lockdown','Post-lockdown','no lockdown','Location','northeast')
xlabel('Time (days)')
ylabel('Number of people')
title('Infection dynamics with lockdown at day 160')
xlim([0 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find final sizes and duration under lockdown

%Duration of outbreak with lockdown
Last_Inft2 = find(Classes4.Ia+Classes4.Is>1,1,'last'); %finds the last remaining infection
Duration2 = ceil(Classes4.t(Last_Inft2))

%Final size without lockdown
Final_size=round(Classes.R(end)+Classes.D(end))

%Final size with lockdown
Final_size2=round(Classes4.R(end)+Classes4.D(end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Run lockdowns varying the time at which they are introduced and store
%final size
for period1=1:250

    mintime=0;
    maxtime=period1;

    %Rerun with no restriction for 1st period (same param), introduce a "mintime"
    [Classes2] = ODE_SIIHR_model(para,ICs,mintime,maxtime);

    %Now define new parameter and ICs for next phase
    param2=struct('beta1',0.3*0.5,'beta2',0.2*0.5,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
    ICs2= struct('S',Classes2.S(end),'Ia',Classes2.Ia(end),'Is',Classes2.Is(end),'H',Classes2.H(end),'D',Classes2.D(end),'R',Classes2.R(end));

    %Define new start and end times for period 2 (with restrictions)
    mintime = period1;
    maxtime = period1+100;

    %Continue simulation following restriction
    [Classes3] = ODE_SIIHR_model(param2,ICs2,mintime,maxtime);

    %Define new initial conditions and parameter for (post restriction) period 3
    param3=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
    ICs3=struct('S',Classes3.S(end),'Ia',Classes3.Ia(end),'Is',Classes3.Is(end),'H',Classes3.H(end),'D',Classes3.D(end),'R',Classes3.R(end));

    %Define start and end times for period 3
    mintime = period1+100;
    maxtime = 10000;
   

    %Continue simulation following restriction ending
    [Classes4] = ODE_SIIHR_model(param3,ICs3,mintime,maxtime);
    
    Final_size3(period1)=round(Classes4.R(end)+Classes4.D(end));

end
[min,days]=min(Final_size3);

%Plot time of introduction against final size
t=[1:250];

figure(3)
plot(t,Final_size3)
xlim([0 250])
xlabel('Time when lockdown starts (days)')
ylabel('Final size of outbreak')
title('Final size of outbreak against time at which lockdown is introduced')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run lockdowns now varying both time at which it is introduced and length of lockdown, storing final size
%%
for period2=30:180 %vary length

    for period1=1:250 %vary introduction day

        mintime=0;
        maxtime=period1;

        %Rerun with no restriction for 1st period (same param), introduce a "mintime"
        [Classes2] = ODE_SIIHR_model(para,ICs,mintime,maxtime);

        %Now define new parameter and ICs for next phase
        param2=struct('beta1',0.3*0.5,'beta2',0.2*0.5,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
        ICs2= struct('S',Classes2.S(end),'Ia',Classes2.Ia(end),'Is',Classes2.Is(end),'H',Classes2.H(end),'D',Classes2.D(end),'R',Classes2.R(end));

        %Define new start and end times for period 2 (with restrictions)
        mintime = period1;
        maxtime = period1+period2;

        %Continue simulation following restriction
        [Classes3] = ODE_SIIHR_model(param2,ICs2,mintime,maxtime);

        %Define new initial conditions and parameter for (post restriction) period 3
        param3=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
        ICs3=struct('S',Classes3.S(end),'Ia',Classes3.Ia(end),'Is',Classes3.Is(end),'H',Classes3.H(end),'D',Classes3.D(end),'R',Classes3.R(end));

        %Define start and end times for period 3
        mintime = period1+period2;
        maxtime = 600;
   

        %Continue simulation following restriction ending
        [Classes4] = ODE_SIIHR_model(param3,ICs3,mintime,maxtime);
    
        %store final size
        Final_size31(period2,period1)=round(Classes4.R(end)+Classes4.D(end));
    end

end

%%
%Collect the variables needed to plot graph and reshape into correct format
Final_size32=Final_size31([30:180],:); 
 t=[1:250];

t2=repelem(t,151);
length=(30:180);
length2=repmat(length,1,250);
fs=reshape(Final_size32,1,[]);
 lockdownlength=[30:180];

 %Plot all three variables using a heat map
figure(4) 
set(0,'defaultaxesfontsize',17)
colormap(jet)
 scatter3(t2,length2,fs,30,fs,'filled');
 ylim([30 180])
 xlabel('Time when lockdown starts')
 ylabel('Length of lockdown')
view(0,90)
 c=colorbar;
 c.Label.String = 'Final size of outbreak';
 c.Label.FontSize=17;
 title('Effects of different lockdowns on the outbreak')