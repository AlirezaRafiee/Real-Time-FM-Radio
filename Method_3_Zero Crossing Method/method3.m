clc
clear
close all
%Zero Crossing Method
fc=89e06;%Center frequency
FrontEndSampleRate=1e06;
f_sample=FrontEndSampleRate;
FrameLength=2^18;
duration=20;
f_desired=50e03; %frequency of audio to play
RTL_Obj= comm.SDRRTLReceiver(...
    'CenterFrequency', fc ,...
    'EnableTunerAGC', true,...
    'SampleRate', FrontEndSampleRate, ...
    'SamplesPerFrame', FrameLength,...
    'OutputDataType', 'double');

filter1=load('filter1_method3.mat','Num');%%%%Filters used in method demodulate FM
filter1=filter1.Num;

filter2=load('filter2_method3.mat','Num');
filter2=filter2.Num;

filter3=load('filter3_method3.mat','Num');
filter3=filter3.Num;

filterdc=load('DCBlock.mat','Num');
filterdc=filterdc.Num;

time=0;
tic
while(time <= duration)
data=step(RTL_Obj);
data=filter(filter1,1,data);%initial Filtering
data=real(data);%To calculate the zero crossing
data=sign(data);
data_new=zeros(1,length(data)+20);%%%%%monostable
for i=1:1:(length(data)-1)%%%%%monostable
    if( data(i) <= 0 && data(i+1) > 0)
            data_new(i)=1;
    elseif(data_new(i) ~= 1)
        data_new(i)=0;
    end
end
data=data_new(1:1:length(data));
Pre=filter(filter3,1,data);% second filtering
Sound=decimate(Pre,round(f_sample/f_desired)); %Decimation for play
Sound=filter(filter3,1,Sound);%third filtering
Sound=filter(filterdc,1,Sound); %%DCBlock maybe use
sound(10*Sound,f_desired);
time=toc;
end
release(RTL_Obj);