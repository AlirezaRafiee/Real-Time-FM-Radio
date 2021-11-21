clc
clear
close all
%Enevelope or Hilbert method
fc=100e06;%Center frequency
FrontEndSampleRate=1e06;
f_sample=FrontEndSampleRate;
FrameLength=2^18;
duration=10;
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

filter3=load('filter3_method1.mat','Num');
filter3=filter3.Num;

filterdc=load('DCBlock.mat','Num');
filterdc=filterdc.Num;

time=0;
tic
while(time <= duration)
data=step(RTL_Obj);
data=filter(filter1,1,data);%initial Filtering
data=diff(data);%differential of data :: Convert to AM
Enevelope=sqrt(data.^2 + (imag(hilbert(real(data)))).^2+(imag(hilbert(imag(data))).^2)); %Calculate Enevelope
Pre=filter(filter2,1,Enevelope);% second filtering
Sound=decimate(Pre,round(f_sample/f_desired));%Decimation for play
Sound=filter(filter3,1,Sound); %third filtering
%Sound=filter(filterdc,1,Sound); %%DCBlock maybe use
sound(4*real(Sound),f_desired);
time=toc;
end
release(RTL_Obj);