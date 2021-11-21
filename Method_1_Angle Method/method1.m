clc
clear
close all
%%Angle Method
fc=89e06;%Center frequency
FrontEndSampleRate=1e06;
f_sample=FrontEndSampleRate;
FrameLength=2^18;
duration=20;
f_desired=50e03;%frequency of audio to play
RTL_Obj= comm.SDRRTLReceiver(...
    'CenterFrequency', fc ,...
    'EnableTunerAGC', true,...
    'SampleRate', FrontEndSampleRate, ...
    'SamplesPerFrame', FrameLength,...
    'OutputDataType', 'double');

filter1=load('filter1_method1.mat');%%%%Filters used in method demodulate FM
filter1=filter1.Num;

filter2=load('filter2_method1.mat');
filter2=filter2.Num;%sos2tf(filter2.SOS);

filter3=load('filter3_method1.mat','Num');
filter3=filter3.Num;

filterdc=load('DCBlock.mat','Num');
filterdc=filterdc.Num;

time=0;
tic
while(time <= duration)
data= step(RTL_Obj);
data=filter(filter1,1,data); %initial Filtering
phase=unwrap(angle(data));%Calculate phase
phase=diff(phase); %Calculate Frequency
phase=filter(filter2,1,phase); % second filtering
Sound=decimate(phase,round(f_sample/f_desired));%Decimation for play
Sound=filter(filter3,1,Sound); %third filtering
%Sound=filter(filterdc,1,Sound);  %%DCBlock maybe use
sound(Sound,f_desired);
time=toc;
end
release(RTL_Obj);