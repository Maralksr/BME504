% Alex Czaja
% Dr. Valero-Cevas, BME504 Neuromuscular Systems
% Implementation of two-fiber type Hill-type model

% Equations of two-element Hill-type model
% According to "Accuracy of gastrocnemius muscles forces in walking and
% running goats predicted by one-element and two-element Hill-type models"
% by Sabrina S.M. Lee et. al. published in Journal of Biomechanics

% Constants borrowed from Ali's code and from above paper
% Making _ be suffix of global variables to be replaced later for Monte
% Carlo simulations
b_ = [0.73, 9.90];          % [beta_slow, beta_fast]
tau_act_ = [34.06, 18.14];  % [tau_act_slow, tau_act_fast]
theta_ = 0;                 % pennation angle, assuming 0;
c_ = 1;                     % for calculating total muscle force from active and passive elements
l_ = 15.9;                  % optimal fascicle length [mm] lateral gastrocnemius from above paper
v_0_ = 2.74;                % maximum shortening velocity from above paper
k_ = 0.29;                  % force-velocity curvature from above paper


% Define activation functions for evaluating fitness
activ_ = 0 : .001 : 1;     % activation ramp of muscle activation from which to get slow and fast fiber activations


% Running model
% Get activations of fibers over ramp activation of muscle
activ_step = activ_(2) - activ_(1);
activ_slow = zeros(1, length(activ_));
activ_fast = zeros(1, length(activ_));
for i = 2 : length(activ_)
    %populate slow and fast fiber activations
    %all activations start at zero, then use a_dot values to proceed
    [a_slow_dot, a_fast_dot] = activation_transfer([activ_(i-1), activ_slow(i-1), activ_fast(i-1)], tau_act_, b_);
    activ_slow(i) = activ_slow(i-1) + a_slow_dot;
    activ_fast(i) = activ_fast(i-1) + a_fast_dot;
end

% Calculate total muscle force generated over activation function


% Total muscle force
% Fm = c(F_f + F_p(l)) * cos(theta)
% Where F_f (F_f-hat in paper) is active component of the muscle fiber force,
% F_p(l) (F_p-hat(l) in paper) is passive component of muscle fiber force,
% constant c reflects the maximum isometric force generated by muscle,
% and theta is pennation angle.
function F_m = total_muscle_force(a_slow, a_fast, l, v_0, k, c, theta, F_applied, direction)
    v = velocity_from_force(F_applied, v_0, k, direction);
    F_f = total_active_force(a_slow, a_fast, l, v, v_0, k);
    F_p = total_passive_force(l);
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
    F_f = (a_slow * force_length_active(l) * force_velocity(v, v_0(1), k(1))) + ...
        (a_fast * force_length_active(l) * force_velocity(v, v_0(2), k(2)));
end

% Transfer functions
% a1_dot + ((1/tau_act1)*(beta1 + (1-beta1)*EMG(t-t_off))) * a1(t) = (1/tau_act1) * EMG(t-t_off)
% a2_dot + ((1/tau_act2)*(beta2 + (1-beta2)*a1(t))) * a2(t) = (1/tau_act2) * a1(t)
% a3_dot + ((1/tau_act3)*(beta3 + (1-beta3)*a3(t))) * a3(t) = (1/tau_act3) * a2(t)
function [a_dot_slow, a_dot_fast] = activation_transfer(activ, tau_act, b)
    % Returns a_dot as vector [a1_dot, a2_dot, a3_dot]
    % where a1_dot is transfer function of the whole muscle, a2_dot is the
    % transfer function of slow fibers, and a3_dot is the transfer function
    % of the fast fibers.
    % Expects tau_act as vector [tau_act1, tau_act2, tau_act3],
    % b represents beta as vector [beta1, beta2, beta3],
    % EMG as scalar value EMG(t-t_off),
    % and a as vector (a1(t), a2(t), a3(t)]
    %a_dot(1) = (1/tau_act(1))*EMG - ((1/tau_act(1))*(b(1)+(1-b(1))*EMG))*a(1);
    a_dot_slow = (1/tau_act(1))*activ(1) - ((1/tau_act(1))*(b(1)+(1-b(1))*activ(2)))*activ(2);
    a_dot_fast = (1/tau_act(2))*activ(2) - ((1/tau_act(2))*(b(3)+(1-b(3))*activ(3)))*activ(3);
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
    assert((direction == 'concentric' || direction == 'eccentric'), ...
        'velocity_from_force: must pass direction as eccentric or concentric');
    if direction == 'concentric'
        v = ((1 - F) * (v_0 * k)) / (k - F);
    else
        v = (2 - 2 * F) / ((1/v_0) + ((3*7.56)/(v_0*k)) - ((2*7.56*F)/(v_0*k)));
    end
end

% Prototype sigmoid function to use in place of fuzzy logic toolbox's
% figmf(x, [a, c]) function
% Current problem: only displays sigmoidal growth between 0 and 1
function y = sigmoid(x, p)
    y = zeros(numel(x));
    a = p(1);
    c = p(2);
    parfor i = 1 : numel(y)
        y(i) = 1 / (1 + exp(-a*(x(i) - c)));
    end
end