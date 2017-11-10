function out = run_corr(x,y);
%function out = run_corr(x,y);
%
% RUN_CORR  --  This function implements the "real-time" running
%               correlation algorithm.
%
% Input Parameters:
%   x  -- the signal to correlate with (i.e., the transmitted signal)
%   y  -- the signal to correlate against (i.e., the received signal)
%
% Output Parameters:
%   out  -- the resulting correlation signal with size equal to
%             length(x) + length(y) - 1

% Written by Mark Bartsch, Winter 2002
% Modification History:
%    8/16/02  --  Added modification history (MB)



% Initialize the output vector to all zeros.  It should have 
%    length(x) + length(y) - 1 elements.
out = zeros(1,length(x)+length(y)-1);

% Make x and y column vectors, regardless of incoming size
x = x(1:end);
y = y(1:length(y));

% Beyond the end of the signal y, the input is zero
y(end+1:end+length(x)) = 0;

% Buffer is just a column vector.  We initialize it to
%    all zeros (with the same size as x) to begin with.
buffer = zeros(length(x),1);

for n = 1:length(out)
    
    % current_sample is the single-sample input to our "correlation box"
    current_sample = y(n);
    
    % Update the buffer.  Discard the first sample in the buffer, and put
    %   current_sample at the end of the buffer.  
    
    buffer(1:end-1) = buffer(2:end,1);
    buffer(end,1) = current_sample;
   
    % From x and the contents of buffer, we can calculate
    %   the next output sample
    
    storage_array = zeros(length(x),1);
    for counter = 1:length(x)
        storage_array(counter,1)=x(counter,1)*buffer(counter,1);
    end
    out(n) = sum(storage_array);

end
