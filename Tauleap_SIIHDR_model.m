%Stochastic (tau leap) SIIHDR model code

function [Classes] = Tauleap_SIIHDR_model(para,ICs,maxtime,tau)

 
%Run the tauleap algorithm for an SIIHDR model

%Store the starting point of the simultion from the ICs and copy to a new
%structure called Classes. Define the starting time to be 0
Classes = ICs;
Classes.t = 0;

%Define the current state of the model from the ICs
S=ICs.S;
Ia=ICs.Ia;
Is=ICs.Is;
H=ICs.H;
R=ICs.R;
D=ICs.D;
%Define the current time in the model as 0
t=0;

%Run the model until either the maxtime is exceeded or until there are no
%infected and infectious people remaining
while ((t<maxtime) && ((Ia>0) || (Is>0)))

%Define event rates 
sym_infection = para.p_S*(para.beta1*Ia+para.beta2*Is)*S/para.N;
asym_infection= (1-para.p_S)*(para.beta1*Ia+para.beta2*Is)*S/para.N;
hospital= para.p_H*para.gamma2*Is;
hosp_recover= (1-para.p_D)*para.xi*H;
death=para.p_D*para.xi*H;
sym_recovery = (1-para.p_H)*para.gamma2*Is;
asym_recovery = para.gamma1*Ia;

%Compute how many events occur for each time step
sym_inf_events=poissrnd(tau*sym_infection);
asym_inf_events=poissrnd(tau*asym_infection);
hospital_events=poissrnd(tau*hospital);
hospital_recovery_events=poissrnd(tau*hosp_recover);
death_events=poissrnd(tau*death);
sym_recovery_events=poissrnd(tau*sym_recovery);
asym_recovery_events=poissrnd(tau*asym_recovery);


%Update events

S = S - sym_inf_events -asym_inf_events;
Is = Is + sym_inf_events - sym_recovery_events - hospital_events;
Ia = Ia + asym_inf_events - asym_recovery_events;
H= H + hospital_events - hospital_recovery_events - death_events;
D = D + death_events;
R = R + sym_recovery_events + asym_recovery_events + hospital_recovery_events;

%Check nothing less than zero
if S<0 %too many people were infected so "undo" infections (both asymptomatic and symptomatic)
        tmp=round(Is/(Ia+Is)*S);
        S=0;
        Is=Is+tmp; %take a proportion back out of Is class
        Ia=Ia+(1-tmp); %take remaining proportion back out of Ia class
end
if Ia<0 %too many people have recovered so "undo" recovery
        tmp=Ia;
        Ia=0;
        R=R+tmp;
end
if Is<0 %too many people have either been hospitalised or recovered "undo" recovery and hospitalisations
        tmp=round(H/(H+R)*Is);
        I=0;
        H=H+tmp; %take a proportion back out of H class
        R=R + (1-tmp);%take remaining back out of recovery class
end
if H<0 %too many people died or recovered
        tmp=round(D/(D+R)*H);
        H=0;
        D=D+tmp; %take a proportion back out of D class
        R=R+(1-tmp); %take remaining proportion back out of R class
end
%Note D and R will always be grater than or equal to 0


%Update time
t = t+tau;

%Save information in the Classes structure by extending the matrices of the
%model state and the associated time
Classes.t = [Classes.t t];
Classes.S = [Classes.S S];
Classes.Ia = [Classes.Ia Ia];
Classes.Is = [Classes.Is Is];
Classes.D = [Classes.D D];
Classes.H = [Classes.H H];
Classes.R = [Classes.R R];
end
    