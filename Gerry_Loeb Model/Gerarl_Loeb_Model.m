L0M= 14.05  %bicep  %optimal fascicle length --> the length at which the muscle produces maximal tetanic force





% replaced with v_0.5 -->velocity of shortening required to produce 0.5F0
% at L0 during a tetanic stimulus -->used to scale the force velocity
% dependencies
%and f_0.5 --> stimulus frequency necessary to produce 0.5 F0 at L0 during
%a tetanic, isometric contraction; used to scale rise and fall time.
%Tau_C % time scaling parameter for maximal muscle shortening velocity and rise and fall times





% calculate maximal force output of the muscle based on PCSA
%alpha = muscle_parameter.pennationAngle; % pennation angle ignore
%mass = muscle_parameter.mass; % muscle mass  %ignore
density = 1.06; % g/cm^3
mass=1
PCSA = (mass*1000)/density/L0M; %physiological cross-sectional area
Specific_Tension = 31.4;
F0M = PCSA * Specific_Tension; % maximal force %Maximal tetanic isometric force
Fmax = 10; % maximal voluntary force determined empilically
Fbaseline =1; % baseline force without muscle activation


%%



%%


%%Active force components
TL= 0.088;
dL_eff=(Muscle_Length-L0M).^3/(TL*(1-A_f));    %effective length- function of L, t, Af


af=0.56; nf0=2.1;  nf1=[5,3.3];
nf=nf0+nf1*(1/L0M-1);
A_f=1-exp(-(Y*S*f_eff/(af*nf))^nf);      %Activation force relationship- function of L_eff, f_eff, S, Y



%%
%slow
Tf1_s=24.2;  Tf2_s=16;  Tf3_s=33.2;  Tf4_s=17.8;
%fast
Tf1_f=20.6;  Tf2_f=13.6;  Tf3_f=28.2;
if df_eff>=0
    Tf= Tf1*L^2+Tf2*f_env;
else
    Tf=(Tf3+Tf4*A_f)/L;
end
df_int=(f_env - f_int)/Tf;   % f_int --> intermediate firing frequency 
df_eff=(f_int - f_eff)/Tf;    %Rise and fall time- function of t, L, Af, f_env--- f_eff --> effective firing frequency

%%


%slow
as1_s=1;  as2_s=1;
%fast
as1_f=1.76;  as2_f=0.96; Ts_f=43;

if f_eff<0.1 
    as=as1;
else
    as=as2;
end
dS= (as-S)/Ts;        %Sag- function of t, f_eff
%%

%Slow
cY_s=0.35;  VY=0.1;  TY=200;
%fast
cY_f=0; 
dY=(1-cY*(1-exp(-abs(V)/VY))-Y)/TY;        %Yield- function of t, V
%%





    


F_CE=F_L *F_V* A_f;   %Active

F_total= F_PE + F_CE;    %total contractile element force