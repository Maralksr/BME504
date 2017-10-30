%% Hill Model

% Model: 4 Element Model. 
% Damper, Activator, and spring (Kp) in
% parallel. Attached to spring (Ks) in series.

% ACTIVATION, A
% FORCE, F
% DAMPING, B
% PARALLEL STIFFNESS, Fp
% SERIES STIFFNES, Fs

close all;

%for Activation
T1 = 20;    %time 'on'
T2 = 14;    %time 'off'
n = 5;      %number of cycles

%for Mechanical Components of Model
%B = 0.4;       %damping constant //de-comment for assuming Damping is
                %constant, and comment out B damping equation
Ks = 0.02;      %series stiffness constant
Kp = 3;         %parallel stiffness constant
Db = 0.005;
Rb = 0.003;

%Activation Equation
A = @(t) ((t/(T1+T2) - floor(t/(T1+T2))) < T1/(T1+T2));
%Damping Constant Equation____ B = @(t) (D_b*(1-R_b)*F/max(A) + R_b)
%solved alongside Force, from Reference 1. Assumption: max(A) = 1

dFBdt = @(t,FB) [ (A(t) - FB(1))*Ks/(FB(2) + Kp);
                 Db*((1-Rb)*FB(1) + Rb)];

%ode inputs for Activation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tspan1 = [0 n*(T1+T2)];

[t, FB] = ode45(dFBdt, tspan1, [0 0]', options);

hold on
plot(t, FB(:,1), 'LineWidth', 3)
plot(t, A(t), 'LineWidth', 2)
plot(t, FB(:,2), 'LineWidth', 2)
title('Output Force for given Activation and Damping Pattern')
legend('Force', 'Activation', 'Damping Pattern')
axis([0 tspan1(2) 0 1.3])
descr = {strcat('t pulse = 16 ms; t interval =', sprintf('%5d', T2), 'ms')};

%% Monte Carlo

%Random number generator for picking values of parameters:
tstart = clock;
seed = round(sum(1000*tstart));
rand('state',seed);             

%iteration numbers
N = 100;

%for Activation
T1 = 20;    %time 'on'
T2 = 14;    %time 'off'
n = 5;      %number of cycles

%initialize Force(1) and Damping (2) by iteration matrix
%FB_matrix = zeros(600,3, N);
time_matrix = zeros(600, N);
Force_matrix = zeros(600, N);
Damp_matrix = zeros(600 ,N);
iterations_matrix = meshgrid(1:100,1:600);

%Activation Equation
A = @(t) ((t/(T1+T2) - floor(t/(T1+T2))) < T1/(T1+T2));

%ode inputs for Activation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tspan1 = [0 n*(T1+T2)];

for i = 1: N

clear t FB

%ranges of parameters
Ks_r = random('unif', 0.01, 0.03);        %series stiffness constant
Kp_r = random('unif', 1, 4);              %parallel stiffness constant
Db_r = random('unif', 0.004, 0.05);       %Damping eq constants
Rb_r = random('unif', 0.003, 0.03);       

%Damping Constant Equation____ B = @(t) (D_b*(1-R_b)*F/max(A) + R_b)
%solved alongside Force, from Reference 1. Assumption: max(A) = 1

dFBdt = @(t,FB) [ (A(t) - FB(1))*Ks_r/(FB(2) + Kp_r);
                 Db_r*((1-Rb_r)*FB(1) + Rb_r)];

[t, FB] = ode45(dFBdt, tspan1, [0 0]', options);

% Build matrix for Time, Force, and Damping
time_matrix(1:length(FB),i) = t;
Force_matrix(1:length(FB),i) = FB(:,1);
Damp_matrix(1:length(FB),i) = FB(:,2);

%FB_matrix(1:length(FB),:,i) = horzcat(t,FB);

end

%plot with for loop
figure()
hold on
mesh(time_matrix,iterations_matrix,Force_matrix)
title('Output Force for given Activation and Damping Pattern')
xlabel('Time')
ylabel('Iteration')
zlabel('Force');

figure()
mesh(time_matrix,iterations_matrix,Damp_matrix)
title('Corresponding Damping Force Pattern')
xlabel('Time')
ylabel('Iteration')
zlabel('Damping Force');






%%
% References

% 1.  Damping Constant Equation
% S. Schmitt, M. Günther, T. Rupp, A. Bayer, and D. Häufle, 
% ?Theoretical Hill-Type Muscle and Stability: Numerical Model and Application,? 
% Computational and Mathematical Methods in Medicine, vol. 2013, 
% Article ID 570878, 7 pages, 2013. doi:10.1155/2013/570878
