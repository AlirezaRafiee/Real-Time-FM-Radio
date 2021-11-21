function varargout = project_matlab(varargin)
% PROJECT_MATLAB MATLAB code for project_matlab.fig
%      PROJECT_MATLAB, by itself, creates a new PROJECT_MATLAB or raises the existing
%      singleton*.
%
%      H = PROJECT_MATLAB returns the handle to a new PROJECT_MATLAB or the handle to
%      the existing singleton*.
%
%      PROJECT_MATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT_MATLAB.M with the given input arguments.
%
%      PROJECT_MATLAB('Property','Value',...) creates a new PROJECT_MATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project_matlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project_matlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project_matlab

% Last Modified by GUIDE v2.5 01-Feb-2021 00:25:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project_matlab_OpeningFcn, ...
                   'gui_OutputFcn',  @project_matlab_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before project_matlab is made visible.
function project_matlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project_matlab (see VARARGIN)

% Choose default command line output for project_matlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes project_matlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project_matlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc

Background=imread('Radio.jpg'); %to set background in start of program
axes(handles.axes2); % Background to axis2 
image(Background); %Show background

Onoff= get(handles.togglebutton2,'Value'); %On / off key

Mute= get(handles.Mute,'Value'); % Mute key
type=get(handles.selecttype,'value'); % Music Player or Radio Type

filter1_1=load('filter1_method1.mat'); %%%%Filters used in method demodulate FM
filter1_1=filter1_1.Num;
filter2_1=load('filter2_method1.mat');
filter2_1=filter2_1.Num;
filter3_1=load('filter3_method1.mat','Num');
filter3_1=filter3_1.Num;
filter1_2=load('filter1_method3.mat','Num');
filter1_2=filter1_2.Num;
filter2_2=load('filter2_method3.mat','Num');
filter2_2=filter2_2.Num;
filter3_2=load('filter3_method1.mat','Num');
filter3_2=filter3_2.Num;
filter1_3=load('filter1_method3.mat','Num');
filter1_3=filter1_3.Num;
filter2_3=load('filter2_method3.mat','Num');
filter2_3=filter2_3.Num;
filter3_3=load('filter3_method3.mat','Num');
filter3_3=filter3_3.Num;
filterdc=load('DCBlock.mat','Num');
filterdc=filterdc.Num;

userdata=load('userInput.mat');%%Data for DSPtoolBox method
userInput=userdata.userInput;

while(Onoff==0) %if Radio is off do nothing
    alert= sprintf('\n The Radio is Off . Turn On it and Try again');
    set(handles.text4, 'string', alert); %show alert to output
    Onoff= get(handles.togglebutton2,'Value');
    drawnow;
end 
if (type == 2) %Music Player Type
    music=get(handles.edit1,'string'); %Music Name
    [music,freq]=audioread(music); %Music Read
    sound(music,freq); %Play music
elseif(type == 1) %Play Radio Type
    
dB_Plot=get(handles.radiobutton4,'Value'); %plot in DB or Not Key

Noise_Reject=get(handles.radiobutton3,'Value');%Use pwletch or Not Key

fc = get(handles.edit1,'string');%Center Frequency 
fc = str2double(fc);

sec = get(handles.edit2,'string'); %Time Duration
sec = str2double(sec);

FMtype=get(handles.FMtype,'value');%Type of method

playplot=get(handles.popupmenu1,'value');%to recognize plot or play or twince
if( playplot == 1 )
    plot_en = 1; %%%if plot_en ==1 the program will show the PSD  
    play =0; %%%If play == 1 the program will play radio
elseif ( playplot == 2 )
    plot_en = 0;
    play =1;
else 
    plot_en = 1 ;
    play = 1 ;
end
FrontEndSampleRate=1e06;%F_sample
f_sample=FrontEndSampleRate; 
FrameLength=2^18;
f_desired=50e03;%frequency of audio to play


RTL_Obj= comm.SDRRTLReceiver(...
    'CenterFrequency', fc ,...
    'EnableTunerAGC', true,...
    'SampleRate', FrontEndSampleRate, ...
    'SamplesPerFrame', FrameLength,...
    'OutputDataType', 'double');

F_axes=linspace(fc-FrontEndSampleRate/2,fc+FrontEndSampleRate/2,FrameLength); %%Frequency axis
time=0;%real time
TR2=0; %maximum ratio  in no noise rejection
CR2=0; %maximum PSD in no noise rejection
TRN=0;%maximum ratio  in noise rejection
CRN=0;%maximum PSD in noise rejection
tic %Start time 
Rec=[0];%For Record Radio

axes(handles.axes1); %%%Choose axes 1 to plot PSD


while( time < sec )
    data= step(RTL_Obj);
    Mute= get(handles.Mute,'Value');
    Record= get(handles.Record,'Value');%%if Record == 1 the program Record radio
if(Noise_Reject==0 && plot_en == 1) %% PSD = (fourier Transform of data )^2
    Data_FFT=fftshift(fft(real(data)))/FrontEndSampleRate;
    PSD=abs(Data_FFT).^2;
    T=round(FrameLength/6); %for FM Recignize
    sumation=sum(PSD((2*T+1):1:(4*T)));
    total=sum(PSD);
    TR=sumation/total;
    CR1=max(PSD);
    CR2=max(CR1,CR2); %find max PSD in whole Duration
    TR2=max(TR,TR2); %find max ratio in whole Duration
    if(dB_Plot==1) % Plot in dB scale
        Data_FFT=10*log(PSD);
        plot(F_axes,Data_FFT,'red');
        xlabel('Frequency(Hz)');
        xlim([F_axes(2*T+1),F_axes(4*T)]);
        ylabel('|Sound(j\Omega)|^2(dB)');
        title('Power Spectrum Density in dB');
        grid on 
        grid minor
    else  %%Plot in normal scale
        plot(F_axes,PSD,'red');
        xlabel('Frequency(Hz)');
        xlim([F_axes(2*T+1),F_axes(4*T)]);
        ylabel('|Sound(j\Omega)|^2');
        title('Power Spectrum Density');
        grid on 
        grid minor
    end
        if (TR2>0.9 || TRN > 0.9 || CR2>= 0.2 || CRN >= 0.2 ) %%For recognize FM
            FM_Rec= sprintf('\n you see a FM Spectrum');
            set(handles.text4, 'String', FM_Rec);   
        else
            FM_Rec= sprintf('\n this center frequency contains no FM Spectrum');
            set(handles.text4, 'String', FM_Rec);
        end
elseif(Noise_Reject==1 && plot_en == 1 ) %%PSD=pwelch
    PSD = pwelch(data,FrameLength,0,F_axes,FrontEndSampleRate);
    T=round(FrameLength/6); %for FM Recignize
    sumation=sum(PSD((2*T+1):1:(4*T)));
    total=sum(PSD);
    TR=sumation/total;
    TRN=max(TR,TRN); %find max ratio in whole Duration
    CR1=max(PSD);
    CRN=max(CR1,CRN); %find max PSD in whole Duration
        if(dB_Plot==1) % Plot in dB scale
        Data_FFT=10*log(PSD);
        plot(F_axes,Data_FFT,'red');
        xlabel('Frequency(Hz)');
        xlim([F_axes(2*T+1),F_axes(4*T)]);
        ylabel('|Sound(j\Omega)|^2(dB)');
        title('Power Spectrum Density in dB By Pwelch');
        grid on 
        grid minor
        else   % Plot in normal scale
        plot(F_axes,PSD,'red');
        xlabel('Frequency(Hz)');
        xlim([F_axes(2*T+1),F_axes(4*T)]);
        ylabel('|Sound(j\Omega)|^2');
        title('Power Spectrum Density By Pwelch');
        grid on 
        grid minor
        end
        if (TR2>0.9 || TRN > 0.9 || CR2>= 0.2 || CRN >= 0.2 ) %%For recognize FM
            FM_Rec= sprintf('\n you see a FM Spectrum');
            set(handles.text4, 'String', FM_Rec);   
        else
            FM_Rec= sprintf('\n this center frequency contains no FM Spectrum');
            set(handles.text4, 'String', FM_Rec);
        end
end
if(FMtype == 1 && play == 1) %%Angle Method
    data=filter(filter1_1,1,data); %initial Filtering
    phase=unwrap(angle(data)); %Calculate phase
    phase=diff(phase); %Calculate Frequency
    phase=filter(filter2_1,1,phase); % second filtering
    Sound=decimate(phase,round(f_sample/f_desired)); %Decimation for play
    Sound=filter(filter3_1,1,Sound);  %third filtering
    %Sound=filter(filterdc,1,Sound);  %%DCBlock maybe use
    if(Mute == 0) %%if mute =1 dont play anything
    sound(Sound,f_desired);
    end
    if(Record == 1) %%if Record ==1 save the Sound 
        Rec=[Rec ; Sound];
    end
elseif(FMtype == 2 && play == 1)%Enevelope or Hilbert method
    data=filter(filter1_2,1,data);  %initial Filtering
    data=diff(data);%differential of data :: Convert to AM
    Enevelope=sqrt(data.^2 + (imag(hilbert(real(data)))).^2+(imag(hilbert(imag(data))).^2)); %Calculate Enevelope
    Pre=filter(filter1_2,1,Enevelope); % second filtering
    Sound=decimate(Pre,round(f_sample/f_desired)); %Decimation for play
    Sound=filter(filter3_2,1,Sound); %third filtering
    Sound=filter(filterdc,1,Sound); %%DCBlock maybe use
    if(Mute == 0) %%if mute =1 dont play anything
    sound(4*real(Sound),f_desired);
    end
    if(Record == 1) %%if Record ==1 save the Sound
        Rec=[Rec ; 4*real(Sound)];
    end
elseif(FMtype == 3 && play == 1) %Zero Crossing Method
    data=filter(filter1_3,1,data); %initial Filtering
    data=real(data); %To calculate the zero crossing
    data=sign(data);
    data_new=zeros(1,length(data)+20); %%%%%monostable
    for i=1:1:(length(data)-1) %%%%%monostable
        if( data(i) <= 0 && data(i+1) > 0)
                data_new(i)=1;
        elseif(data_new(i) ~= 1)
            data_new(i)=0;
        end
    end 
    data=data_new(1:1:length(data));
    Pre=filter(filter3_3,1,data); % second filtering
    Sound=decimate(Pre,round(f_sample/f_desired));  %Decimation for play
    Sound=filter(filter3_3,1,Sound); %third filtering
    %Sound=filter(filterdc,1,Sound); %%DCBlock maybe use
    if(Mute == 0) %%if mute =1 dont play anything
    sound(10*Sound,f_desired);
    end
    if(Record == 1)  %%if Record ==1 save the Sound
        Rec=[Rec ; Sound];
    end
elseif(FMtype == 4 && play == 1) %%DSPtoolBox method
    release(RTL_Obj); %%Release RTL :: because we want to set a new one
    userInput.CenterFrequency=fc;% Set CenterFrequency
    % Calculate FM system parameters based on the user input
    [fmRxParams,sigSrc] = helperFMConfig(userInput);

    % Create FM broadcast receiver object and configure based on user input
    fmBroadcastDemod = comm.FMBroadcastDemodulator(...
    'SampleRate', fmRxParams.FrontEndSampleRate, ...
    'FrequencyDeviation', fmRxParams.FrequencyDeviation, ...
    'FilterTimeConstant', fmRxParams.FilterTimeConstant, ...
    'AudioSampleRate', fmRxParams.AudioSampleRate, ...
    'Stereo', false);
    %   Create audio player
    player = audioDeviceWriter('SampleRate',fmRxParams.AudioSampleRate);
    % Initialize radio time
    radioTime = time;
    % Main loop
    while radioTime < sec
        Record= get(handles.Record,'Value');
        % Receive baseband samples (Signal Source)
        if fmRxParams.isSourceRadio
            if fmRxParams.isSourcePlutoSDR
            rcv = sigSrc();
            lost = 0;
            late = 1;
        else
            [rcv,~,lost,late] = sigSrc();
            end
        else
            rcv = sigSrc();
            lost = 0;
            late = 1;
        end
        % Demodulate FM broadcast signals and play the decoded audio
         audioSig = fmBroadcastDemod(rcv);
         if(Mute == 0) %%if mute =1 dont play anything
         player(audioSig);
         end
        if(Record == 1) %%if Record ==1 save the Sound
        Rec=[Rec ; double(audioSig) ];
        end
        % Update radio time. If there were lost samples, add those too.
        radioTime = radioTime + fmRxParams.FrontEndFrameTime + ...
        double(lost)/fmRxParams.FrontEndSampleRate;
        time_taken=sprintf('time : %1.1d (sec)',radioTime); % for Show real time in output
        set(handles.text3, 'String',time_taken);
        drawnow;
    end
end
    time_taken=sprintf('time : %1.1d (sec)',time); % for Show real time in output
    set(handles.text3, 'String',time_taken);
    drawnow;
    time=toc;
end
audiowrite('RecordFM.WAV',Rec,f_desired); % save the Recored audio in 'RecordFM.WAV'
end
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton2.
function radiobutton2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
Background=imread('Radio.jpg'); % to set and Show the Background
axes(handles.axes2);
image(Background);

% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Mute.
function Mute_Callback(hObject, eventdata, handles)
% hObject    handle to Mute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mute


% --- Executes on selection change in selecttype.
function selecttype_Callback(hObject, eventdata, handles)
% hObject    handle to selecttype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selecttype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selecttype


% --- Executes during object creation, after setting all properties.
function selecttype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selecttype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FMtype.
function FMtype_Callback(hObject, eventdata, handles)
% hObject    handle to FMtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FMtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FMtype


% --- Executes during object creation, after setting all properties.
function FMtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FMtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Record.
function Record_Callback(hObject, eventdata, handles)
% hObject    handle to Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Channels.
function Channels_Callback(hObject, eventdata, handles)
% hObject    handle to Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Channels


% --- Executes during object creation, after setting all properties.
function Channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_channel.
function Save_channel_Callback(hObject, eventdata, handles)
% hObject    handle to Save_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Onoff= get(handles.togglebutton2,'Value'); %On / off key
while(Onoff==0) %if Radio is off do nothing
    alert= sprintf('\n The Radio is Off . Turn On it and Try again');
    set(handles.text4, 'string', alert); %show alert to output
    Onoff= get(handles.togglebutton2,'Value');
    drawnow;
end 

fc = get(handles.edit1,'string');
fc = str2double(fc); %center frequency we want to save
channel_number=get(handles.Channels,'value'); %to find which of the channels we choose for save
Channel_memory=load('Channels.mat','Channels');%read memory
Channel_amount=Channel_memory.Channels;
x(1)=Channel_amount(1);
x(2)=Channel_amount(2);
x(3)=Channel_amount(3);
x(4)=Channel_amount(4);
x(5)=Channel_amount(5);
x(channel_number)=fc; %write in specific memory location
Channels=x;
save('Channels.mat','Channels'); % write the memory


% --- Executes on button press in Read_Channel.
function Read_Channel_Callback(hObject, eventdata, handles)
% hObject    handle to Read_Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Onoff= get(handles.togglebutton2,'Value'); %On / off key
while(Onoff==0) %if Radio is off do nothing
    alert= sprintf('\n The Radio is Off . Turn On it and Try again');
    set(handles.text4, 'string', alert); %show alert to output
    Onoff= get(handles.togglebutton2,'Value');
    drawnow;
end 

channel_number=get(handles.Channels,'value'); %to find which of the channels we choose for read
Channel_memory=load('Channels.mat'); %read memory
Channel_amount=Channel_memory.Channels; %read memory
x(1)=Channel_amount(1); 
x(2)=Channel_amount(2);
x(3)=Channel_amount(3);
x(4)=Channel_amount(4);
x(5)=Channel_amount(5);
fc=x(channel_number);
fc_write=sprintf('%d',fc); %Set the center frequency 
set(handles.edit1,'String',fc_write); %show the center frequency in 'center frequency block'
drawnow;


% --- Executes on button press in Searcher_forward.
function Searcher_forward_Callback(hObject, eventdata, handles)
% hObject    handle to Searcher_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Onoff= get(handles.togglebutton2,'Value'); %On / off key
while(Onoff==0) %if Radio is off do nothing
    alert= sprintf('\n The Radio is Off . Turn On it and Try again');
    set(handles.text4, 'string', alert); %show alert to output
    Onoff= get(handles.togglebutton2,'Value');
    drawnow;
end 

fc = get(handles.edit1,'string');
fc = str2double(fc); %center frequency we want to save
TR=0;
CR1=0;
FrontEndSampleRate=1e06;%F_sample
FrameLength=2^18;
k=0;
while( TR < 0.8 &&  CR1 < 0.2 && fc < 108e06 )
    k=k+1;
    fc=fc+100e03;
    RTL_Obj= comm.SDRRTLReceiver(...
    'CenterFrequency', fc ,...
    'EnableTunerAGC', true,...
    'SampleRate', FrontEndSampleRate, ...
    'SamplesPerFrame', FrameLength,...
    'OutputDataType', 'double');
    data= step(RTL_Obj);
    Data_FFT=fftshift(fft(real(data)))/FrontEndSampleRate;
    PSD=abs(Data_FFT).^2;
    T=round(FrameLength/6); %for FM Recignize
    sumation=sum(PSD((2*T+1):1:(4*T)));
    total=sum(PSD);
    TR=sumation/total;
    CR1=max(PSD); 
    message= sprintf('\n Searching Forward %d -th try',k);
    set(handles.text4, 'String', message); 
    drawnow;
    release(RTL_Obj);
end
message= sprintf('\n Search Finished');
set(handles.text4, 'String', message);
fc_write=sprintf('%d',fc); %Set the center frequency 
set(handles.edit1,'String',fc_write); %show the center frequency in 'center frequency block'
drawnow;


% --- Executes on button press in Searcher_Backward.
function Searcher_Backward_Callback(hObject, eventdata, handles)
% hObject    handle to Searcher_Backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Onoff= get(handles.togglebutton2,'Value'); %On / off key
while(Onoff==0) %if Radio is off do nothing
    alert= sprintf('\n The Radio is Off . Turn On it and Try again');
    set(handles.text4, 'string', alert); %show alert to output
    Onoff= get(handles.togglebutton2,'Value');
    drawnow;
end 

fc = get(handles.edit1,'string');
fc = str2double(fc); %center frequency we want to save
TR=0;
CR1=0;
FrontEndSampleRate=1e06;%F_sample
FrameLength=2^18;
k=0;
while(TR < 0.8 &&  CR1 < 0.2 && fc > 88e06 )
    fc=fc-100e03;
    RTL_Obj= comm.SDRRTLReceiver(...
    'CenterFrequency', fc ,...
    'EnableTunerAGC', true,...
    'SampleRate', FrontEndSampleRate, ...
    'SamplesPerFrame', FrameLength,...
    'OutputDataType', 'double');
    data= step(RTL_Obj);
    Data_FFT=fftshift(fft(real(data)))/FrontEndSampleRate;
    PSD=abs(Data_FFT).^2;
    T=round(FrameLength/6); %for FM Recignize
    sumation=sum(PSD((2*T+1):1:(4*T)));
    total=sum(PSD);
    TR=sumation/total;
    CR1=max(PSD); 
    message= sprintf('\n Searching Backward %d -th try',k);
    set(handles.text4, 'String', message); 
    drawnow;
    release(RTL_Obj);
end
message= sprintf('\n Search Finished');
set(handles.text4, 'String', message);
fc_write=sprintf('%d',fc); %Set the center frequency 
set(handles.edit1,'String',fc_write); %show the center frequency in 'center frequency block'
drawnow;
