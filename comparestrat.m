%Code used to run comparisons in section 4


%%Clear all the previous variables and close previous figures
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First run outbreak with vaccinations at day 90
%Define model parameters as a structure
parav = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0,'alpha',0.61,'ro',0); 

%Define initial conditions as a structure
ICs = struct('S',parav.N-1,'Sv',0,'Ia',0,'Is',1,'H',0,'D',0,'R',0);

%Define time to run model for
mintime = 0;
period1=90;
maxtime = period1;

%Run model by calling function ODE_SIIHRv_model.m
[ClassesV] = ODE_SIIHRv_model(parav,ICs,mintime,maxtime);

% Now define new parameter and ICs for next phase
paramv2=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0.005,'alpha',0.61,'ro',0.005); 
ICs2= struct('S',ClassesV.S(end),'Sv',ClassesV.Sv(end),'Ia',ClassesV.Ia(end),'Is',ClassesV.Is(end),'H',ClassesV.H(end),'D',ClassesV.D(end),'R',ClassesV.R(end));

%Define new start and end times for period 2 (with restrictions)
mintime = period1;
maxtime = 800;

%Continue simulation following restriction
[ClassesV2] = ODE_SIIHRv_model(paramv2,ICs2,mintime,maxtime);

[Classes3] = ODE_SIIHR_model(parav,ICs,0,maxtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lockdown
%Define start and end times
periodL1=167;
mintimeL=0;
maxtimeL=periodL1;

paraL = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);

%Rerun with no restriction for 1st period (same param), introduce a "mintime"
[ClassesL2] = ODE_SIIHR_model(paraL,ICs,mintimeL,maxtimeL);

%Now define new parameter and ICs for next phase
paramL2=struct('beta1',0.3*0.5,'beta2',0.2*0.5,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
ICsL2= struct('S',ClassesL2.S(end),'Ia',ClassesL2.Ia(end),'Is',ClassesL2.Is(end),'H',ClassesL2.H(end),'D',ClassesL2.D(end),'R',ClassesL2.R(end));

%Define new start and end times for period 2 (with restrictions)
mintimeL = periodL1;
maxtimeL = periodL1+120;

%Continue simulation following restriction
[ClassesL3] = ODE_SIIHR_model(paramL2,ICsL2,mintimeL,maxtimeL);

%Define new initial conditions and parameter for (post restriction) period 3
paramL3=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000);
ICsL3=struct('S',ClassesL3.S(end),'Ia',ClassesL3.Ia(end),'Is',ClassesL3.Is(end),'H',ClassesL3.H(end),'D',ClassesL3.D(end),'R',ClassesL3.R(end));

%Define start and end times for period 3
mintimeL = periodL1+120;
maxtimeL = 10000;

%Continue simulation following restriction ending
[ClassesL4] = ODE_SIIHR_model(paramL3,ICsL3,mintimeL,maxtimeL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quarantine
%Define model parameters as a structure
paraQ = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'d_1',0.05,'tau',0.1); 

%Define initial conditions as a structure
ICsQ = struct('S',paraQ.N-1,'Ia',0,'Is',1,'H',0,'D',0,'R',0,'Q',0);

%Define time to run model for
mintimeq = 0;
maxtimeq = 1000;

%Run model by calling function ODE_SIIHRmodel.m
[ClassesQ] = ODE_SIIHRq_model(paraQ,ICsQ,mintimeq,maxtimeq);

%Plot the dynamics
figure(2)
set(0,'defaultaxesfontsize',24)
set(0,'defaultlinelinewidth',3)

set(0,'defaultaxesfontsize',16)

plot(ClassesV.t,ClassesV.Ia+ClassesV.Is,'r',ClassesV2.t,ClassesV2.Ia+ClassesV2.Is,'r')
hold on
plot(ClassesL2.t,ClassesL2.Ia+ClassesL2.Is,'b',ClassesL3.t,ClassesL3.Ia+ClassesL3.Is,'b',ClassesL4.t,ClassesL4.Ia+ClassesL4.Is,'b')
hold on
plot(ClassesQ.t,ClassesQ.Ia+ClassesQ.Is,'g')

legend('Vaccinations','','Lockdown','','','Quarantine','','northeast')
xlabel('Time (days)')
ylabel('Number of people')
xlim([0 500])
ylim([0 inf])
title('Infection dynamics under different interventions')

%Find final size and duration for all 3 interventions

%Duration of outbreak with lockdown
Last_InftL = find(ClassesL4.Ia+ClassesL4.Is>1,1,'last'); %finds the last remaining infection
DurationL = ceil(ClassesL4.t(Last_InftL));

%Duration of outbreak with vaccine
Last_InftV = find(ClassesV2.Ia+ClassesV2.Is>1,1,'last'); %finds the last remaining infection
DurationV = ceil(ClassesV2.t(Last_InftV));

%Duration of outbreak with quarantine
Last_InftQ = find(ClassesQ.Ia+ClassesQ.Is>1,1,'last'); %finds the last remaining infection
DurationQ = ceil(ClassesQ.t(Last_InftQ));

%Final size with lockdown
deathsL=round(ClassesL4.D(end))

%Final size with vaccine
deathsV=round(ClassesV2.D(end))

deathsQ=round(ClassesQ.D(end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finding metrics for CEA

%Set up time horizon
NumYears = 10;
maxtime = 365*NumYears + 1;

%Gives start point (index) of each year
t_Yr=[1:365:NumYears*365+1];


%Set up the number of strategies
NumStrats = 4;

%Preallocate matrices for storing DALYs and Costs
DALYMat = zeros(1,NumStrats);

% Run the models for each startegy for the number of iterations and store the key
% metrics

for Strat = 1:NumStrats
   
    if Strat==1
        %run no intervention for 10 years
        maxtime=365*10;
        [Classes3] = ODE_SIIHR_model(parav,ICs,0,maxtime);

               %Compute person-years and deaths each year
        for r=1:length(t_Yr)-1
            %Person years infected
            PersonYears_annual(r) = trapz(Classes3.t([t_Yr(r): t_Yr(r+1)]), Classes3.Is([t_Yr(r): t_Yr(r+1)]))/365; %divide by 365 to go from person days to person years
    
            %Deaths
             Deaths_annual(r) = Classes3.D(t_Yr(r+1))- Classes3.D(t_Yr(r));
  
        end
    
        %Compute DALYs and store
        DALYs = dw*PersonYears_annual + yll*Deaths_annual;
      
       
    elseif Strat==2
        parav = struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0,'alpha',0.61,'ro',0); 

        %Define initial conditions as a structure
        ICs = struct('S',parav.N-1,'Sv',0,'Ia',0,'Is',1,'H',0,'D',0,'R',0);

        %Define time to run model for
        mintime = 0;
        period1=100;
        maxtime = period1;

        %Run model by calling function ODE_SIIHRv_model.m
        [ClassesV] = ODE_SIIHRv_model(parav,ICs,mintime,maxtime);

        % Now define new parameter and ICs for next phase
        paramv2=struct('beta1',0.3,'beta2',0.2,'p_S',0.75,'gamma1',1/5,'gamma2',1/10,'p_H',0.338,'xi',0.125,'p_D',0.235,'N',68000000,'u',0.005,'alpha',0.61,'ro',0.005); 
        ICs2= struct('S',ClassesV.S(end),'Sv',ClassesV.Sv(end),'Ia',ClassesV.Ia(end),'Is',ClassesV.Is(end),'H',ClassesV.H(end),'D',ClassesV.D(end),'R',ClassesV.R(end));

        %Define new start and end times for period 2 (with restrictions)
        mintime = period1;
        maxtime = 365*10;

        %Continue simulation following restriction
        [ClassesV2] = ODE_SIIHRv_model(paramv2,ICs2,mintime,maxtime);
        
        combineV_t = [ClassesV.t(1:end-1);ClassesV2.t];
        combineV_Is= [ClassesV.Is(1:end-1);ClassesV2.Is];
        combineV_D = [ClassesV.D(1:end-1);ClassesV2.D];
        
                       %Compute person-years and deaths each year
        for r=1:length(t_Yr)-1
            %Person years infected
            PersonYears_annual(r) = trapz(combineV_t([t_Yr(r): t_Yr(r+1)]), combineV_Is([t_Yr(r): t_Yr(r+1)]))/365; %divide by 365 to go from person days to person years
    
            %Deaths
             Deaths_annual(r) = combineV_D(t_Yr(r+1))- combineV_D(t_Yr(r));
  
        end
    
        %Compute DALYs and store
        DALYs = dw*PersonYears_annual + yll*Deaths_annual;

    elseif Strat==3
    %Define start and end times
    periodL1=167;
    mintimeL=0;
    maxtimeL=periodL1;

    %Rerun with no restriction for 1st period (same param), introduce a "mintime"
    [ClassesL2] = ODE_SIIHR_model(paraL,ICs,mintimeL,maxtimeL);

    %Continue simulation following restriction
    [ClassesL3] = ODE_SIIHR_model(paramL2,ICsL2,mintimeL,maxtimeL);

    %Define start and end times for period 3
    mintimeL = periodL1+120;
    maxtimeL = 365*10;

    %Continue simulation following restriction ending
    [ClassesL4] = ODE_SIIHR_model(paramL3,ICsL3,mintimeL,maxtimeL);

    combineL_t=[ClassesL2.t(1:end-1);ClassesL3.t(1:end-1);ClassesL4.t];
    combineL_Is=[ClassesL2.Is(1:end-1);ClassesL3.Is(1:end-1);ClassesL4.Is];
    combineL_D=[ClassesL2.D(1:end-1);ClassesL3.D(1:end-1);ClassesL4.D];

                              %Compute person-years and deaths each year
        for r=1:length(t_Yr)-1
            %Person years infected
            PersonYears_annual(r) = trapz(combineL_t([t_Yr(r): t_Yr(r+1)]), combineL_Is([t_Yr(r): t_Yr(r+1)]))/365; %divide by 365 to go from person days to person years
    
            %Deaths
             Deaths_annual(r) = combineL_D(t_Yr(r+1))- combineL_D(t_Yr(r));
  
        end
    
        %Compute DALYs and store
        DALYs = dw*PersonYears_annual + yll*Deaths_annual;

    elseif Strat==4
    %Define time to run model for
    mintimeq = 0;
    maxtimeq = 365*10;

    %Run model by calling function ODE_SIIHRq_model.m
    [ClassesQ] = ODE_SIIHRq_model(paraQ,ICsQ,mintimeq,maxtimeq);

                          %Compute person-years and deaths each year
        for r=1:length(t_Yr)-1
            %Person years infected
            PersonYears_annual(r) = trapz(ClassesQ.t([t_Yr(r): t_Yr(r+1)]), ClassesQ.Is([t_Yr(r): t_Yr(r+1)]))/365; %divide by 365 to go from person days to person years
    
            %Deaths
             Deaths_annual(r) = ClassesQ.D(t_Yr(r+1))- ClassesQ.D(t_Yr(r));
  
        end
    
        %Compute DALYs and store
        DALYs = dw*PersonYears_annual + yll*Deaths_annual;
    end

        
        DALYMat(1,Strat) = sum(DALYs,2);
       
        
        DALYs_annual{Strat}(1,:)=DALYs;
end 
  %%
%Calculate delta DALYs  
DDALYMat = repmat(DALYMat(:,1),1,NumStrats) - DALYMat
DCost =[0, 11.7e9, 83.3e9, 110.4e9]
%ICERs show lockdown is dominated - recalculate ICERs excluding
%lockdown
ICER2=zeros(1,3);
ICER2(1)=0;
ICER2(2)=(DCost(2)-DCost(1))/(DDALYMat(2)-DDALYMat(1));
ICER2(3)=(DCost(4)-DCost(2))/(DDALYMat(4)-DDALYMat(2));

%New mean delta cost vector for table
NewDCost=zeros(1,3);
NewDCost(1)=DCost(1);
NewDCost(2)=DCost(2);
NewDCost(3)=DCost(4);
%New  delta DALY vector 
NewDDALY=zeros(1,3);
NewDDALY(1)=DDALYMat(1);
NewDDALY(2)=DDALYMat(2);
NewDDALY(3)=DDALYMat(4);

%New updated ICER table
ICERtable2 = array2table([NewDCost' NewDDALY' ICER2'])
ICERtable2.Properties.VariableNames={'Delta Costs','Delta DALYs','ICER'};





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot CE plane
figure(3)
clf
set(0,'defaultaxesfontsize',18)
set(0,'defaultlinelinewidth',2)
hold on
%Draw on ICERs
ax = gca; 
ax.ColorOrderIndex = 2; 
plot([0 NewDDALY(2)], [0 NewDCost(2)]./1e9,'--')
hold on
plot([NewDDALY(2) NewDDALY(3)], [NewDCost(2) NewDCost(3)]./1e9,'--')

ax = gca; 
ax.ColorOrderIndex = 1; %resets colour cycle to number 2 (red)
h(1)=scatter(0,0,100,'o','filled');
h(2)=scatter(NewDDALY(2),NewDCost(2)./1e9,100,'o','filled');
h(3)=scatter(NewDDALY(3),NewDCost(3)./1e9,100,'o','filled');
axis([0 Inf 0 Inf])
ylabel('Additional costs (billions £)')
xlabel('DALYs averted')
title('CE plane for 3 strategies over 10 year time horizon')
legend(h,'No intervention','Vaccinations','Quarantine','location','NorthWest')