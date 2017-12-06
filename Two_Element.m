%% Alex Czaja
% Dr. Valero-Cevas, BME504 Neuromuscular Systems
% Implementation of two-fiber type Hill-type model

% Equations of two-element Hill-type model
% According to "Accuracy of gastrocnemius muscles forces in walking and
% running goats predicted by one-element and two-element Hill-type models"
% by Sabrina S.M. Lee et. al. published in Journal of Biomechanics

format compact;
close all;
clear all;

%% Global/Initial Parameters
% Constants borrowed from Ali's code and from above paper
% Making _ be suffix of global variables to be replaced later for Monte
% Carlo simulations

% [beta_slow, beta_fast]
b_ = [0.73, 9.90];

% [tau_act_slow, tau_act_fast]
tau_act_ = [34.06, 18.14];

% pennation angle
theta_ = 0;

% for calculating total muscle force from active and passive elements
% consider renaming to fatigue factor becuase c
% should conceptually be 0 <= c <= 1
c_ = 1;

% optimal fascicle length lateral gastrocnemius from above paper
% given in mm, converted to m
l_opt_ = 15.9e-2;

% maximum shortening velocity from above paper
v_0_ = 2.74;

% [k_slow, k_fast] force-velocity curvature from above paper
k_ = [0.18, 0.29];

% Set ode solver tolerance params
%options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);


%% Running Basic Model
% Define activation functions for evaluating fitness
%dt = 0.001;
%t = 0 : dt : 10;
t_ = 0 : 1000;  % [ms]
activ_ = sigmoid(t_, [0.025, 500]);

[t, a, ~] = ode45(@(t, a) dadt(t, a, tau_act_, b_), [0, 1000], [0, 0]');
activ_slow = a(:, 1)';
activ_fast = a(:, 2)';

% Calculate total muscle force generated over activation function
force = zeros(1, size(a, 1));
for i = 1 : size(a, 1)
    %total_muscle_force(a_slow, a_fast, l, l_opt, v_0, k, c, theta, F_applied, direction)
    %muscle_force(a_slow, a_fast, l, l_opt, v, v_0, k, c, theta)
    %force(i) = total_muscle_force(activ_slow(i), activ_fast(i), l_opt_, l_opt_, v_0_, k_, c_, theta_, 0, 'isometric');
    force(i) = muscle_force(activ_slow(i), activ_fast(i), l_opt_, l_opt_, 0, v_0_, k_, c_, theta_);
end

% Plot a and force
figure;
plot(t_, activ_, 'b');
hold on;
plot(t, activ_slow, 'g');
hold on;
plot(t, activ_fast, 'r');
title('Activation Curves of Fibers');
legend('Whole Muscle Activation', 'Slow Fiber Activation', 'Fast Fiber Activation', ...
    'Location', 'Southeast');
ylabel('Activation Level');
xlabel('Time [ms]');
legend('Whole muscle', 'Slow fibers', 'Fast fibers');

figure;
plot(t, force);
title('force');
title('Force Production');
ylabel('Force [N]');
xlabel('Time [ms]');

%% Monte Carlo
% For two-element model, want to stochastically model
% beta1, beta2, tau1, tau2, k1, and k2

% Initialize rng
tstart = clock;
seed = round(sum(1000*tstart));
rand('state', seed);

% Define maximum number of iterations
MAX_ITER = 100;

% Preallocate matrices to store results, keep structs of params used for
% each iter, and keep vectors of force produced

% M x N: M -> structs of params, N -> which iteration struct was made
MC_params = cell(1, MAX_ITER);

% M x N x P: M -> [a_slow, a_fast], N -> a(t), P -> which iteration
MC_activ = zeros(2, 1, MAX_ITER);

% M x N: M -> which iteration, N -> times corresponding to soln from solver
MC_time = cell(1, MAX_ITER);

% M x N: M -> which iteration, N -> force(t)
MC_force = zeros(MAX_ITER, 1);

% Perform Monte Carlo
for i = 1 : MAX_ITER
    % Generate params and store in MC_params
    params.b = [
        random('unif', 0.5*b_(1), 1.5*b_(1)), ...
        random('unif', 0.5*b_(2), 1.5*b_(2))];
    params.tau = [
        random('unif', 0.5*tau_act_(1), 1.5*tau_act_(1)), ...
        random('unif', 0.5*tau_act_(2), 1.5*tau_act_(2))];
    params.k = [
        random('unif', 0.5*k_(1), 1.5*k_(1)), ...
        random('unif', 0.5*k_(2), 1.5*k_(2))];
    MC_params{i} = params;
    
    % Integrate activation for fibers using new params
    [t, a] = ode45(@(t, a) dadt(t, a, params.tau, params.b), [0, 1000], [0, 0]');
    a = a';
    t = t';
    
    % Store solution
    % Recall MC_activ: M x N x P: M -> [a_slow, a_fast], N -> a(t), P -> which iteration
    % Recall MC_time: M x N: M -> which iteration, N -> times corresponding to soln from solver
    MC_activ(1:size(a, 1), 1:size(a, 2), i) = a;
    MC_time{i} = t;
    
    % Use activation soln to calculate force(t)
    force = zeros(1, length(t));
    for j = 1 : length(force)
        force(j) = muscle_force(a(1, j), a(2, j), l_opt_, l_opt_, 0, v_0_, params.k, c_, theta_);
    end
    
    % Store solution
    % Recall MC_force: M x N: M -> which iteration, N -> force(t)
    MC_force(i, 1:length(force)) = force;
end

% SANITY CHECK
% Check what fraction of MC trials produced activation curves where fast
% fibers activate more quickly than slow fibers. This would be inaccurate
% behavior.
% for i = 1 : MAX_ITER
%     figure;
%     plot(MC_time{i}, MC_activ(1, 1:length(MC_time{i}), i), 'g');
%     hold on;
%     plot(MC_time{i}, MC_activ(2, 1:length(MC_time{i}), i), 'r');
% end
% title('MC activation');


% Find where the maximum force occurred after MC modelling
% Recall form above that M dim of MC_force corresponds to iteration number
% max(mat, [], 2) returns a column vector of the max of each row of mat
% and max(vec) returns [max_of_vec, idx_of_max], where idx_of_max is the
% row corresponding to the MC iteration that produced the best results
[f, which_iter] = max(max(MC_force, [], 2));
fprintf("Max force produced by MC: %d\n", f);
fprintf("Iter where MC produced max force: %i\n", which_iter);

% Plot force production after MC modelling contained in MC_force
figure;
subplot(1, 2, 1);
for i = 1 : MAX_ITER
    plot(MC_time{i}, MC_force(i, 1:length(MC_time{i})));
    hold on;
end
title('MC Force(t)');
ylabel('Force [N]');
xlabel('Time');

% Include maximal force curve in subplot
subplot(1, 2, 2);
plot(MC_time{which_iter}, MC_force(which_iter, 1:length(MC_time{which_iter})));
title('Max Force(t)');
ylabel('Force [N]');
xlabel('Time');

% Pull best parameter set and activation curve the produced maximum force
best_params = MC_params{which_iter};
best_activ = MC_activ(:, :, which_iter);

%% Concentric and Eccentric Modelling

% Weight of mass attached to muscle [N]
%applied_forces = [0.50, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 10.0];
%applied_forces = [1.0, 2.0, 3.0];
applied_forces = [0.5];

% Weightings of activation curve to test for each weight
%attenuations = [0.2, 0.4, 0.6, 0.8, 0.95];
attenuations = [.33, .66, .98];

% Perform modelling for each applied weight for each activation attenuation
% Create simulation to show in a figure for a single applied weight with a
% curve for simulation using each attenuation
%[t, y] = ode15s(@(t, y) dYdt(t, y, best_params, l_opt_, v_0_, c_, theta_, 1, applied_force), [1, 1000], [0, 0]');
for app = 1 : length(applied_forces)
    fprintf('Modelling force = %.2f N\n', applied_forces(app));
    % Store displacement, velocity, time the model ran for, and fiber 
    % activations that were used during that time interval for calculating
    % force produced
    modeled_v = cell(1, length(attenuations));
    modeled_y = cell(1, length(attenuations));
    modeled_t = cell(1, length(attenuations));
    modeled_a = cell(1, length(attenuations));
    modeled_f = cell(1, length(attenuations));
    
    % Calculate displacement and plot
    parfor att = 1 : length(attenuations)
        % Copy broadcast variables for parfor loop
        params_copy = best_params;
        applied_forces_copy = applied_forces;
        
        fprintf('\tWith a = %.2f\n', attenuations(att));
        % Determine a_slow and a_fast for this attenuation attempt for this
        % weight to pass to dYdt for integration
        dt = 0.1;
        temp_t = 0 : dt : 1000;
        sig = sigmoid(temp_t, [0.025, 500]);
        sig(sig > attenuations(att)) = 0;
        [~, idx] = max(sig);
        att_t = idx * dt;
        
        fprintf('\t\tGetting fiber activations\n');
        [~, a] = ode45(@(t_a, a) dadt(t_a, a, params_copy.tau, params_copy.b), [0, att_t], [0, 0]');
        a_slow = a(end, 1);
        a_fast = a(end, 2);
        modeled_a{att} = [a_slow, a_fast];
        init_force = muscle_force(a_slow, a_fast, l_opt_, l_opt_, 0, v_0_, k_, c_, theta_);
        
        fprintf('\t\tIntegrating F = ma\n');
        % Integrate the solution using attenuated activations
        [t, y] = ode45(@(t, y) dYdt(t, y, a_slow, a_fast, params_copy, l_opt_, v_0_, c_, theta_, applied_forces_copy(app)), [0, 10], [0, 0]', options);
        
        % TEST: use euler method instead of ode solver and save result
%         dt = 1 / 100000;
%         t = [0 : dt : 5];
%         y = zeros(length(t), 2);
%         for j = 1 : length(t)-1
%             dy = dYdt(t(j), y(j, :), a_slow, a_fast, params_copy, l_opt_, v_0_, c_, theta_, applied_forces_copy(app));
%             y(j+1, :) = y(j, :) + (dt * dy');
%         end
        
        modeled_y{att} = y(:, 1);
        modeled_v{att} = y(:, 2);
        %modeled_f{att} = y(:, 3);
        modeled_t{att} = t;
        
    end
    
    figure;
    subplot(1, 2, 1);
    for m = 1 : length(attenuations)
        hold on;
        plot(modeled_t{m}, modeled_y{m});
    end
    title(sprintf('Modeled Movement When\n%.1f N Weight Applied', applied_forces(app)));
    ylabel('Displacement [cm]');
    xlabel('Time');
    %legend('Att=0.2', 'Att=0.4', 'Att=0.6', 'Att=0.8', 'Att=1.0', 'Location', 'Southeast');
    legend('a=0.3', 'a=0.6', 'a=0.95', 'Location', 'Southwest');
    
    for m = 1 : length(attenuations)
        f = zeros(1, length(modeled_t{m}));
        for j = 1 : size(f, 2)
            f(j) = muscle_force(modeled_a{m}(1), modeled_a{m}(2), modeled_y{m}(j), l_opt_, modeled_v{m}(j), v_0_, k_, c_, theta_);
        end
        modeled_f{m} = f;
    end
    subplot(1, 2, 2);
    for i = 1 : length(attenuations)
        hold on;
        plot(modeled_t{i}, modeled_f{i});
    end
    title(sprintf('Modeled Force when\n%.1f N Weight Applied', applied_forces(app)));
    ylabel('Force');
    xlabel('Time');
    %legend('Att=0.2', 'Att=0.4', 'Att=0.6', 'Att=0.8', 'Att=1.0', 'Location', 'Southeast');
    legend('a=0.3', 'a=0.6', 'a=0.95', 'Location', 'Southwest');
end

% Insert inverted modelling and plotting strategy from deprecated_functions
% if need to produce inverse figures of multiple forces per activation

%% Model Functions

% Function for integrating motion instead of current anon fxn dydt
function [out] = dYdt(t, y, a_slow, a_fast, params, l_opt, v_0, c, theta, applied_force)
    %disp(t);
    %https://www.mathworks.com/matlabcentral/answers/92701-how-do-i-pass-out-extra-parameters-using-ode23-or-ode45-from-the-matlab-ode-suite
    f = muscle_force(a_slow, a_fast, l_opt+y(1), l_opt, y(2)/l_opt, v_0, params.k, c, theta);
    out = [ ...
        y(2)*l_opt; ...
        %(applied_force - muscle_force(a_slow, a_fast, l_opt+y(1), l_opt, y(2)/l_opt, v_0, params.k, c, theta)) / applied_force/9.8];
        (applied_force - f) / (applied_force/9.8)]; ...
        %f];
end

% Function for integrating activation
function a = dadt(t, y, tau, b)
    a = zeros(2, 1);
    a(1) = (1/tau(1))*sigmoid(t, [0.025, 500]) - ((1/tau(1))*(b(1)+(1-b(1))*sigmoid(t, [0.025, 500])))*y(1);
    a(2) = (1/tau(2))*y(1) - ((1/tau(2))*(b(2)+(1-b(2))*y(1)))*y(2);
end

% Total muscle force that doesn't depend on passing contraction condition
% and depends on v explicitly. Better for integration during simulation.
function F_m = muscle_force(a_slow, a_fast, len, l_opt, v, v_0, k, c, theta)
    F_f = total_active_force(a_slow, a_fast, len/l_opt, v, v_0, k);
    F_p = force_length_passive(len/l_opt);
    F_m = c * (F_f + F_p) * cos(theta);    
end

% Active component of muscle force determined by this two-element model
% Force is the sum of the force from the slow fibers and the force from the
% fast fibers
% Equations for individual component of force is
% F_n = a(t) * F_a(l) * F_v(v)
function F_f = total_active_force(a_slow, a_fast, l, v, v_0, k)
    % Expects a as vector of activation levels [a_slow(t), a_fast(t)],
    % l as current fascicle length, and v as current fiber velocity,
    % v_0 as vector [v_0_slow, v_0_fast], and k as vector [k_slow, k_fast]
    %l = le3;
    %disp(force_length_active(l));
    F_f = (a_slow * force_length_active(l) * force_velocity(v, v_0, k(1))) + ...
        (a_fast * force_length_active(l) * force_velocity(v, v_0, k(2)));
end

% Force-length relationship (F_a-hat(l) and F_p-hat(l) in paper)
% F_a(l) = (-878.24*(l * 1.253)^2 + 2200.4*(l * 1.254) - 1192) / 186.24
function F_active = force_length_active(len)
    % Returns force-length relationship F_l of active component
    % l is fascicle length
    F_active = (-878.25*(len*1.253)^2 + 2200.4*(len*1.254) - 1192) / 186.24;
end

% Force-length relationship (F_p-hat(l) in paper)
% F_p(l) = exp(-1.3 + 3.8*(l * 1.253)) / 186.24
function F_passive = force_length_passive(len)
    % Returns force-length relationship F_l of passive component
    % l is fascicle length
    F_passive = exp(-1.3 + 3.8*(len*1.253)) / 186.24;
end

% Force-velocity relationship (F_v-hat(v) in paper)
% F_v(v) = { (1 - (v/v_0)) / (1 + (v/(v_0*k)))                  if v <= 0
%          { 1.5 - 0.5*((1 + (v/v_0)) / (1 - ((7.5*v)/(v_0*k))) if v > 0
function F_v = force_velocity(v, v_0, k)
    % Returns force-velocity relationship F_v as scalar F_v-hat
    % Expects v as v(t), v_0 as maximum unloaded shortening velocity,
    % and k as force-velocity curvature
    v = v / v_0;
    if v <= 0 
        F_v = (1 - (v/v_0)) / (1 + (v/(v_0*k)));
    else
        F_v = 1.5 - 0.5 * ((1 + (v/v_0)) / (1 - ((7.56*v)/(v_0*k))));
    end
end

%% Helper Functions

% Prototype sigmoid function to use in place of fuzzy logic toolbox's
% figmf(x, [a, c]) function
% For use later to simulate fast increase in activation level
function y = sigmoid(x, p)
    y = 1 ./ (1 + exp(-p(1) .* (x - p(2))));
end