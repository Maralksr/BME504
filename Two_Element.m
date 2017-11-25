%% Alex Czaja
% Dr. Valero-Cevas, BME504 Neuromuscular Systems
% Implementation of two-fiber type Hill-type model

% Equations of two-element Hill-type model
% According to "Accuracy of gastrocnemius muscles forces in walking and
% running goats predicted by one-element and two-element Hill-type models"
% by Sabrina S.M. Lee et. al. published in Journal of Biomechanics

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
l_opt_ = 15.9e-3;

% maximum shortening velocity from above paper
v_0_ = 2.74;

% [k_slow, k_fast] force-velocity curvature from above paper
k_ = [0.18, 0.29];

% Set ode solver tolerance params
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

%% Running model
% Define activation functions for evaluating fitness
%dt = 0.001;
%t = 0 : dt : 10;
t_ = 0 : 1000;  % [ms]
activ_ = sigmoid(t_, [0.025, 500]);

% Activation transfer function as an anon fxn for ode integration
dadt = @(t, y, tau, b) [ ...
    (1/tau(1))*sigmoid(t, [0.025, 500]) - ((1/tau(1))*(b(1)+(1-b(1))*sigmoid(t, [0.025, 500])))*y(1); ...
    (1/tau(2))*y(1) - ((1/tau(2))*(b(2)+(1-b(2))*y(1)))*y(2)];

% Get activations of fibers over ramp activation of muscle
% activ_step = activ_(2) - activ_(1);
% activ_slow = zeros(1, length(activ_));
% activ_fast = zeros(1, length(activ_));
% for i = 2 : length(activ_)
%     %populate slow and fast fiber activations
%     %all activations start at zero, then use a_dot values to proceed to populate a(t) for fiber types
%     [a_slow_dot, a_fast_dot] = activation_transfer([activ_(i-1), activ_slow(i-1), activ_fast(i-1)], tau_act_, b_);
%     activ_slow(i) = activ_slow(i-1) + a_slow_dot;
%     activ_fast(i) = activ_fast(i-1) + a_fast_dot;
% end

[t, a] = ode45(@(t, a) dadt(t, a, tau_act_, b_), [0, 1000], [0, 0]');
activ_slow = a(:, 1)';
activ_fast = a(:, 2)';

% Calculate total muscle force generated over activation function
force = zeros(1, size(a, 1));
for i = 1 : size(a, 1)
    force(i) = total_muscle_force(activ_slow(i), activ_fast(i), l_opt_, l_opt_, v_0_, k_, c_, theta_, 0, 'isometric');
end

% TEMP: plot a and force
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

%SANITY CHECKS
%MC_force_slow = zeros(1, length(t));
%MC_force_fast = zeros(1, length(t));

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
        force(j) = total_muscle_force( ...
            a(1, j), a(2, j), l_opt_, l_opt_, v_0_, params.k, c_, 0, theta_, 'isometric');
    end
    
    % Store solution
    % Recall MC_force: M x N: M -> which iteration, N -> force(t)
    MC_force(i, 1:length(force)) = force;
end

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
title('MC force');

% Include maximal force curve in subplot
subplot(1, 2, 2);
plot(MC_time{which_iter}, MC_force(which_iter, 1:length(MC_time{which_iter})));
title('Max Force');

% Pull best parameter set and activation curve the produced maximum force
best_params = MC_params{which_iter};
best_activ = MC_activ(:, :, which_iter);

%% Concentric Modelling

% Copy activation and force of best iter
activation = MC_activ(:, :, which_iter);
force = MC_force(which_iter, :);
a_slow = activation(1,:);
a_fast = activation(2,:);

% Weight of mass attached to muscle [N]
applied_force = .25;
%applied_weights = [0.25, 0.50, 0.75, 1, 1.25, 1.5, 2, 2.25];

% Weightings of activation curve to test for each weight
%activation_weights = [0.25, 0.50, 0.75, 1];

% Set up force balance scenario of second order diff eq as system of first
% order diff eqs. (y_d -> y_dot; y_dd -> y_dot_dot)
dydt = @(t, y) [ ...
    y(1); ...
    %total_active_force(a_slow(t), a_fast(t), l, y(2), v_0, k)];
    %(total_active_force(.7, .3, l_opt_, y(2), v_0_, best_params.k) - applied_force) / applied_force/9.8];
    % Recall total_active_force(a_slow, a_fast, l, v, v_0, k)
    (total_active_force(a_slow(end), a_fast(end), ...
        l_opt_-y(2), y(1), v_0_, best_params.k) - applied_force) / applied_force/9.8];
    % function F_m = total_muscle_force_diff(a_slow, a_fast, l, l_opt, v, v_0, k, c, theta, F_applied, direction)
    %(total_muscle_force_diff(activation(1, end), activation(2, end), ...
    %    l_opt_+(sum(y(:,1).*sum(t))), l_opt_, y(1), v_0_, best_params.k, c_, theta_, applied_force, 'concentric') - applied_force) / (applied_force/9.8)];

[t, y] = ode15s(dydt, [0, 10], [0, 0]', options);
figure;
plot(t, y);
title('concentric');

%% Eccentric Modelling


%% Model Functions

% Total muscle force
% Fm = c(F_f + F_p(l)) * cos(theta)
% Where F_f (F_f-hat in paper) is active component of the muscle fiber force,
% F_p(l) (F_p-hat(l) in paper) is passive component of muscle fiber force,
% constant c reflects the maximum isometric force generated by muscle,
% and theta is pennation angle.
function F_m = total_muscle_force(a_slow, a_fast, l, l_opt, v_0, k, c, theta, F_applied, direction)
    v = velocity_from_force(F_applied, v_0, k, direction);
    F_f = total_active_force(a_slow, a_fast, l/l_opt, v, v_0, k);
    F_p = force_length_passive(l);
    F_m = c * (F_f + F_p) * cos(theta);
end


% Total muscle force for ode solver where you can pass current velocity
% from previous iter of solver where y(1) is y_dot
function F_m = total_muscle_force_diff(a_slow, a_fast, l, l_opt, v, v_0, k, c, theta, F_applied, direction)
    %v = velocity_from_force(F_applied, v_0, k, direction);
    F_f = total_active_force(a_slow, a_fast, l/l_opt, v, v_0, k);
    F_p = force_length_passive(l);
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

function F_f = total_active_force_slow(a_slow, l, v, v_0, k)
    F_f = a_slow * force_length_active(l) * force_velocity(v, v_0, k(1));
end

function F_f = total_active_force_fast(a_fast, l, v, v_0, k)
    F_f = a_fast * force_length_active(l) * force_velocity(v, v_0, k(2));
end

% Transfer functions
% a1_dot + ((1/tau_act1)*(beta1 + (1-beta1)*EMG(t-t_off))) * a1(t) = (1/tau_act1) * EMG(t-t_off)
% a2_dot + ((1/tau_act2)*(beta2 + (1-beta2)*a1(t))) * a2(t) = (1/tau_act2) * a1(t)
% a3_dot + ((1/tau_act3)*(beta3 + (1-beta3)*a3(t))) * a3(t) = (1/tau_act3) * a2(t)
function [a_dot_slow, a_dot_fast] = activation_transfer(a, tau, b)
    % Returns a_dot as vector [a2_dot, a3_dot]
    % where a1_dot is transfer function of the whole muscle, a2_dot is the
    % transfer function of slow fibers, and a3_dot is the transfer function
    % of the fast fibers.
    % Expects tau_act as vector [tau_act_slow, tau_act_fast],
    % b represents beta as vector [beta_slow, beta_fast],
    % and a as vector (a_whole(t), a_slow(t), a_fast(t)]
    %a_dot(1) = (1/tau_act(1))*EMG - ((1/tau_act(1))*(b(1)+(1-b(1))*EMG))*a(1);
    %not including a_dot of whole muscle in order to use arbitrary
    %activation pattern for testing
    a_dot_slow = (1/tau(1))*a(1) - ((1/tau(1))*(b(1)+(1-b(1))*a(1)))*a(2);
    a_dot_fast = (1/tau(2))*a(2) - ((1/tau(2))*(b(2)+(1-b(2))*a(2)))*a(3);
end

% Force-length relationship (F_a-hat(l) and F_p-hat(l) in paper)
% F_a(l) = (-878.24*(l * 1.253)^2 + 2200.4*(l * 1.254) - 1192) / 186.24
function F_active = force_length_active(l)
    % Returns force-length relationship F_l of active component
    % l is fascicle length
    F_active = (-878.25*(l*1.253)^2 + 2200.4*(l*1.254) - 1192) / 186.24;
end

% Force-length relationship (F_p-hat(l) in paper)
% F_p(l) = exp(-1.3 + 3.8*(l * 1.253)) / 186.24
function F_passive = force_length_passive(l)
    % Returns force-length relationship F_l of passive component
    % l is fascicle length
    F_passive = exp(-1.3 + 3.8*(l*1.253)) / 186.24;
end

% Force-velocity relationship (F_v-hat(v) in paper)
% F_v(v) = { (1 - (v/v_0)) / (1 + (v/(v_0*k)))                  if v <= 0
%          { 1.5 - 0.5*((1 + (v/v_0)) / (1 - ((7.5*v)/(v_0*k))) if v > 0
function F_v = force_velocity(v, v_0, k)
    % Returns force-velocity relationship F_v as scalar F_v-hat
    % Expects v as v(t), v_0 as maximum unloaded shortening velocity,
    % and k as force-velocity curvature
    if v <= 0 
        F_v = (1 - (v/v_0)) / (1 + (v/(v_0*k)));
    else
        F_v = 1.5 - 0.5 * ((1 + (v/v_0)) / (1 - ((7.56*v)/(v_0*k))));
    end
end

% Force from velocity calculated by inverting force-velocity relationship
% above for the purpose of determine new length
% Requires knowledge of concentric or eccentric muscle activity because
% negative velocity corresponds to concentric or shortening movement
function v = velocity_from_force(F, v_0, k, direction)
    % For isometric case, just pass F as 0 because direction should be
    % passed as 'isometric' yielding desired output of velocity = 0 since
    % muscle is not contracting or elongating
    assert((strcmp(direction, 'concentric') || ...
        strcmp(direction, 'eccentric') || ...
        strcmp(direction, 'isometric')), ...
        'velocity_from_force: must pass direction as eccentric or concentric');
    % Using mean curvature of slow and fast fibers because not specified
    % how to use them separately in the papers researched.
    if strcmp(direction, 'concentric')
        v = ((1 - F) * (v_0 * mean(k))) / (mean(k) - F);
    elseif strcmp(direction, 'eccentric')
        v = (2 - 2 * F) / ((1/v_0) + ((3*7.56)/(v_0*mean(k))) - ((2*7.56*F)/(v_0*mean(k))));
    else
        v = 0;
    end
end

%% Helper Functions

% Prototype sigmoid function to use in place of fuzzy logic toolbox's
% figmf(x, [a, c]) function
% For use later to simulate fast increase in activation level
function y = sigmoid(x, p)
    %len = length(x);
    %y = zeros(1, len);
    %a = p(1);
    %c = p(2);
    y = 1 ./ (1 + exp(-p(1) .* (x - p(2))));
    %for i = 1 : len
    %    y(i) = 1 / (1 + exp(-a * (x(i)-c)));
    %end
end