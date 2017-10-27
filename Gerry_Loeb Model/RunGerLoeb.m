%%------------
%Maral Kasiri
%Reference: Cheng J. et al 2000
%Virtual model
%%-----
clear all;
close all;
clc;


%%
% Parameters
%--------------------
% Architectural parameters of muscle of interest (TA in this study) (Elias et al. 2014;Arnold et al. 2010)
Input.pennationAngle = 0; % pennation angle (radians)
Input.mass = 0.15; % muscle mass (grams)
Input.optimalLength = 6.8; % muscle optimal length (cm)
Input.tendonLength = 24.1; % tendon length (cm)
Input.muscleInitialLength = 6.8; % initial muscle length (cm)
Input.tendonInitialLength = 24.1; % initial tendon length (cm)
Input.MVIC = 363; % maximal voluntary isometric force
Input.baseline = 2.043; % passive force at zero activation

%%
% Force Trajectory
%-----------------------
Fs = 10000;
t = 0:1/Fs:10; % Simulation time
%T_MVC = 0.2; %Target MVC
%Trajectory = [zeros(1,1*Fs) 1/2*t(1:2*Fs) ones(1,length(t)-3*Fs)]; % target force trajectory
Ueff=[zeros(1,floor(length(t)/2)) 0.5*ones(1,floor(length(t)/2)+1)];  %activation

%%
%Run the muscle model 
%-----------------------
%See the GerLoeb() function for details
ModelOutput= GerLoeb(Input,Ueff);


% plot output
figure()
plot(t,ModelOutput.Force)

%hold on
%plot(t,Ftarget)
%legend('Output Force','Target Force')
xlabel('Time(sec)','Fontsize',14)
ylabel('Force(N)','Fontsize',14)