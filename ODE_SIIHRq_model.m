%ODE SIIHDR model with quarantine

function [Classes] = ODE_SIIHRq_model(para,ICs,mintime,maxtime)


%Run ODE using ODE45

%Change ODE options to increase RelTol so that numerical and anlaystical
%results match (to nearest person)
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIIHRqmodel, [mintime:1: maxtime], [ICs.S ICs.Ia ICs.Is ICs.H ICs.D ICs.R ICs.Q], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'Ia',pop(:,2),'Is',pop(:,3),'H',pop(:,4),'D',pop(:,5),'R',pop(:,6),'Q',pop(:,7),'t',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIIHRqmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
Ia=pop(2);
Is=pop(3);
H=pop(4);
D=pop(5);
R=pop(6);
Q=pop(7);

%Write down the ODE system
dS = -(para.beta1*S*Ia/para.N + para.beta2*S*Is/para.N);
dIa= (1-para.p_S)*(para.beta1*S*Ia/para.N + para.beta2*S*Is/para.N) - para.gamma1*Ia;
dIs= para.p_S*(para.beta1*S*Ia/para.N + para.beta2*S*Is/para.N) - para.gamma2*Is-para.d_1*Is;
dH= para.p_H*para.gamma2*Is - para.xi*H;
dD= para.p_D*para.xi*H;
dR = para.gamma1*Ia + (1-para.p_H)*para.gamma2*Is +(1-para.p_D)*para.xi*H + para.tau*Q;
dQ=para.d_1*Is - para.tau*Q;

%Reshape the derivatives into a column vector
dPop = [dS; dIa; dIs; dH; dD; dR; dQ];

end

end