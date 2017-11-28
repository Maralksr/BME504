
close all;
clear all;

tau_act_ = [34.06, 18.14];
b_ = [0.73, 9.90];
attenuation = 0.3;

dt = 0.1;
t = 0 : dt : 1000;
sig = attenuation.*sigmoid(t, [.025, 500]);
figure;
plot(t, sig);

% sig(sig > attenuation) = 0;
% figure;
% plot(t, sig);

[~, idx] = max(sig);
att_t = idx * dt

[t, a] = ode45(@(t, a) dadt(t, a, tau_act_, b_), [0, att_t], [0, 0]');
figure;
plot(t, a);





function a = dadt(t, y, tau, b)
    a = zeros(2, 1);
    a(1) = (1/tau(1))*sigmoid(t, [0.025, 500]) - ((1/tau(1))*(b(1)+(1-b(1))*sigmoid(t, [0.025, 500])))*y(1);
    a(2) = (1/tau(2))*y(1) - ((1/tau(2))*(b(2)+(1-b(2))*y(1)))*y(2);
end 

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

function y = sigmoid(x, p)
    y = 1 ./ (1 + exp(-p(1) .* (x - p(2))));
end