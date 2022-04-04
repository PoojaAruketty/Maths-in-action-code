%ODE SI model code

function [Classes] = ODE_SIS_model(para,ICs,maxtime)


%Run ODE using ODE45

%Change ODE options to increase RelTol so that numerical and anlaystical
%results match (to nearest person)
opts = odeset('RelTol',1e-5);
%Only give outputs at defined timepoints
[t, pop] = ode45(@diff_SIRmodel, [0:1: maxtime], [ICs.S ICs.I], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'I',pop(:,2),'t',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
I=pop(2);


%Write down the ODE system
dS = para.gamma*I-para.beta*S*I/para.N ;
dI = para.beta*S*I/para.N - para.gamma*I;

%Reshape the derivatives into a column vector
dPop = [dS; dI];

end

end
