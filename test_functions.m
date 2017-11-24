
close all;

tau = [34.06, 18.14];
b = [0.73, 9.90];

x = 0 : .01 : 10;
y = sigmoid(x, [2, 4]);
plot(x, y);

% Convert activation_transfer from below into anon fxn that can be
% integrated easily by ode numerical solver
dadt = @(t, y) [ ...
    (1/tau(1))*sigmoid(t, [2, 4]) - ((1/tau(1))*(b(1)+(1-b(1))*sigmoid(t, [2, 4])))*y(1); ...
    (1/tau(2))*y(1) - ((1/tau(2))*(b(2)+(1-b(2))*y(1)))*y(2)];
[t, y] = ode45(dadt, [0, 1000], [0, 0]');
plot(t, y);


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
    % In place of Matlab's fuzzy logic toolbox sigmf(x, [a, c]) function
    len = length(x);
    y = zeros(1, len);
    a = p(1);
    c = p(2);
    parfor i = 1 : len
        y(i) = 1 / (1 + exp(-a * (x(i)-c)));
    end
end