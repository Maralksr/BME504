function [QQ] = BasicFilter(data,Fs,Hz,filter_order,type)
%type can be 'highpass', 'lowpass', 'notch', 'bandpass'
%NOTE: if bandpass, syntax for Hz is [lowestHz  highestHz]
%also, order is usually between 2 and 8.  Usually 4 is good.
%rows are different channels
%columns are time
QQ=zeros(size(data));

for i=1:size(data,1)
    
a=[];b=[];filt1=[];
Nqst=Fs/2;

if strcmp(type,'highpass')
[b,a] = butter(filter_order,Hz/Nqst,'high');%
QQ(i,:) = filtfilt(b,a,data(i,:));

elseif strcmp(type,'lowpass')
[b,a] = butter(filter_order,Hz/Nqst,'low');%
QQ(i,:) = filtfilt(b,a,data(i,:));

elseif strcmp(type,'bandpass')
    %first do a low pass up to the highest freq in band
[b,a] = butter(filter_order,Hz(2)/Nqst,'low');%
filt1 = filtfilt(b,a,data(i,:)); 
a=[];b=[];
%then do a high pass starting from the lowest freq in band
[b,a] = butter(filter_order,Hz(1)/Nqst,'high');%
QQ(i,:) = filtfilt(b,a,filt1);

elseif strcmp(type,'notch')
Q= 25;%this usually works for 60 hz
[b, a]=iirnotch(Hz/Nqst,Hz/Nqst/Q);
QQ(i,:) = filtfilt(b,a,data(i,:));
end


end