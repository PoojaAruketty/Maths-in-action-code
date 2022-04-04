%ODE SIIHDR model with vaccination code

function [Classes] = ODE_SIIHRv_model(para,ICs,mintime,maxtime)


%Run ODE using ODE45

%Change ODE options to increase RelTol so that numerical and anlaystical
%results match (to nearest person)
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIIHRvmodel, [mintime:1: maxtime], [ICs.S ICs.Sv ICs.Ia ICs.Is ICs.H ICs.D ICs.R], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'Sv',pop(:,2),'Ia',pop(:,3),'Is',pop(:,4),'H',pop(:,5),'D',pop(:,6),'R',pop(:,7),'t',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIIHRvmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
Sv=pop(2);
Ia=pop(3);
Is=pop(4);
H=pop(5);
D=pop(6);
R=pop(7);

%Write down the ODE system
dS = -(para.beta1*S*Ia/para.N + para.beta2*S*Is/para.N) -para.u*S;
dSv= para.u*S - (1-para.alpha)*(para.beta1*Sv*Ia/para.N + para.beta2*Sv*Is/para.N) - para.ro*Sv;
dIa= (1-para.p_S)*(para.beta1*S*Ia/para.N + para.beta2*S*Is/para.N)+ (1-para.p_S)*(1-para.alpha)*(para.beta1*Sv*Ia/para.N + para.beta2*Sv*Is/para.N) - para.gamma1*Ia;
dIs= para.p_S*(para.beta1*S*Ia/para.N + para.beta2*S*Is/para.N) +para.p_S*(1-para.alpha)*(para.beta1*Sv*Ia/para.N + para.beta2*Sv*Is/para.N) - para.gamma2*Is;
dH= para.p_H*para.gamma2*Is - para.xi*H;
dD= para.p_D*para.xi*H;
dR = para.ro*Sv + para.gamma1*Ia + (1-para.p_H)*para.gamma2*Is +(1-para.p_D)*para.xi*H;

%Reshape the derivatives into a column vector
dPop = [dS; dSv; dIa; dIs; dH; dD; dR];

end

end