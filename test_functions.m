
close all;

x = 0 : .01 : 10;
y = sigmoid(x, [2, 4]);
plot(x, y);

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