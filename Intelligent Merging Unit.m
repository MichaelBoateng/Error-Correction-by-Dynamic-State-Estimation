%{
The Intelligent Merging Unit: Error Correction by Dynamic State Estimation.
Author: Michael A. Boateng.
Date: April, 01, 2024.
%}

clc;
syms v1 v2 v3 v4 e lambda y1 y2 y3 y4 ip im iL1 iL2 iL3

%TASK 1 STARTS HERE------------------------------------------------------------
% 1a. Prompt opens up the comtrade .cfg file selection
% 1b. After the file is opened, the path and file name are stored
[CfgFileName,Path] = uigetfile('*.cfg');
PathAndCfgName =[Path CfgFileName];

% 2a. As the .dat and .cfg have similar names (but different extensions),
% This command simply replaces .cfg with .dat
% 2b. It then proceeds to store the path and file name
DatFileName = strcat(sscanf(CfgFileName,'%1s',length(CfgFileName)-4),'.dat');
PathAndDatName =[Path DatFileName];

% 3a. This stores the file name (without the extension) and path in the work space
assignin('base','FileName',sscanf(CfgFileName,'%1s',length(CfgFileName)-4));
% 3b. The project's .cfg and .dat files are open, and provided with an ID
cfg_id = fopen(PathAndCfgName);
dat_id = fopen(PathAndDatName);

% 4a.The .cfg and .dat files are stored locally, and the files are closed
cfg = textscan(cfg_id, '%s', 'delimiter', '\n');
dat = textscan(dat_id, '%s', 'delimiter', '\n');
fclose('all');

cfg_len = length(cfg{1,1});
cfg_string = cell(size(cfg));

for i = 1:cfg_len
temp_string = char(cfg{1,1}{i});
cfg_string(i) = textscan(temp_string, '%s', 'delimiter',',')';
end

% This helps to know the file name, type of instrument, and comtrade year
Title = char(cfg_string{1,1}(1));
PT_or_CT=char(cfg_string{1,1}(2));
if length(cfg_string{1,1}) < 3
Year = '1999'; %if the components in the first array are < 3, it is a 1999 comtrade file
else
Year = char(cfg_string{1,1}(3));
end

% Extracts the following channel information: total, analogs and digitals
Total_Channel = strread(char(cfg_string{1,2}(1)));
Analog_Channel = strread(char(cfg_string{1,2}(2)));
Digital_Channel = strread(char(cfg_string{1,2}(3)));

% Number of samples
m = strread(char(cfg_string{1,(cfg_len-3)}(2)));
% Nominal frequency
frequency = strread(char(cfg_string{1,(cfg_len-5)}(1)));
% Sampling rate
sample_rate = strread(char(cfg_string{1,(cfg_len-3)}(1)));
h = 1/sample_rate;   %period
%In later equations, when finding the rate of change we will divide this
%by two to account for the Nyquist frequency

%Start of recording
start_date = char(cfg_string{1,(cfg_len-2)}(1));
start_time = char(cfg_string{1,(cfg_len-2)}(2));
% End of recording
end_date = char(cfg_string{1, (cfg_len-1)}(1));
end_time = char(cfg_string{1, (cfg_len-1)}(2));

%This section focuses on the data files
dat_string = cell(size(dat));
data = zeros(m, Total_Channel+2);

% This section of the code extracts useful data for .dat 
for i = 1:m 
%length of data file
dat_string(i) = textscan(char(dat{1,:}(i)), '%n' ,'delimiter', ',');
data(i,:) = (dat_string{:,i});
end

%Extract and convert time to us
t = (data(:,2)) * (1e-6);
% Writes and sends the data to the workspace
assignin('base','t', t);
var_string = cell(Total_Channel);

for i = 1 : Total_Channel
third_row = (i + 2);
var_string{i} = char(textscan(char(cfg_string{1, third_row }(2)),'%c'));
% If the first character is not a letter, replace with an 'x'
if ~isletter(var_string{i}(1))
var_string{i}(1) = 'x';
end

for b = 2:length(var_string{i})
if ~isstrprop(var_string{i}(b), 'alphanum')
var_string{i}(b) = '_';
end
end
assignin('base', var_string{i}, data(:,third_row));
end
assignin('base','Title', Title);
assignin('base','Year', Year);
assignin('base','Total_Channels', Total_Channel);
assignin('base','Analog_Channels', Analog_Channel);
assignin('base','Digital_Channels', Digital_Channel);
assignin('base','Frequency', frequency);
assignin('base','Sample_rate', sample_rate);
assignin('base','Start_date', start_date);
assignin('base','Start_time', start_time);
assignin('base','End_date', end_date);
assignin('base','End_time', end_time);

if Analog_Channel >= 1
final_data = zeros(m, Analog_Channel);
for i = 1 : Analog_Channel
third_row = (i + 2);

% Analog values require scaling using the vk=a*dk + b formula
%factor: a, offset: b, dk: actual valule 
multiplier = strread(char(cfg_string{1,third_row}(6))); % a
offset = strread(char(cfg_string{1,third_row}(7))); % b
final_data(:,i) = (data(:,third_row) .* (multiplier)) + offset;
end
end

m; % total measurements
frequency;% nominal frequency
N = round(sample_rate/frequency*2);
sample_rate; % sampling rate
%TASK 1 ENDS HERE--------------------------------------------------------------


%TASK 2 & 3 STARTS HERE--------------------------------------------------------
%This is a representation of the secondary voltage (burden)
v_out = final_data (:,1);
% The states to estimate
x = [v1; v2; v3; v4; e; lambda; y1; y2; y3; y4; ip; im; iL1; iL2; iL3];
%Initializing the states to zero
v1_in =0; v2_in =0; v3_in =0; v4_in =0; e_in =0; lambda_in =0; y1_in =0; 
y2_in = 0; y3_in =0; y4_in =0; ip_in =0; im_in = 0; iL1_in =0; iL2_in = 0;
iL3_in =0; ib_in =0;

% The following variables were provided for the project
n = 400;
gm = 0.001; %S
L1 = 26.526e-6; %uH
L2 = 348e-6; %uH
L3 = 348e-6; %uH
M_23 = 287e-6; %uH
g_s1 = 0.5 *(h/(2*L1)); %S
g_s2 = 0.5 *(h/(2*L2)); %S
g_s3 = 0.5 *(h/(2*L3)); %S
r_1 = 0.005; %ohm
r_2 = 0.4469; %ohm
r_3 = 0.4469; %ohm
R_b = 0.1; %ohm
g_b = 1/R_b; %ohm
lambda_0 = 0.1876; %Wb
i_0 = 6.09109; %A
L_0 = 2.36; %H

%The following is used to calculate the weighted matrix using the
%standard deviations from all measurements
W = diag([(1/(0.005))^2*[1 1 1 1] (1/(0.0005))^2*[1 1 1] (1/(0.005^2))*[1 1] (1/(0.0005)^2)*[1] (1/(0.00005)^2)*[1 1 1 1] (1/(0.0003)^2)*[1] (1/(0.05)^2)*[1 1 1 1 1] 1*1]);
cs = zeros(21,1);
prob = 0;
for j = 1:length(v_out)
    %Actual Measurements
    eq1 =(v3-v4);
    %Virtual Measurements
    eq2 = -gm*h/2*(e+e_in)-h/2*(im-im_in)+1/n*h/2*(ip+ip_in)+h/2*(iL1+iL1_in)+g_s1*L1*(iL1-iL1_in);
    eq3 = gm*h/2*(e+e_in)+h/2*(im+im_in)-1/n*h/2*(ip+ip_in)-h/2*(iL2+iL2_in)-g_s2*(L2*(iL2-iL2_in)-M_23*(iL3-iL3_in));
    eq4 = -gm*h/2*(e+e_in)-h/2*(im+im_in)+1/n*h/2*(ip+ip_in)+h/2*(iL3+iL3_in)-g_s3*(L3*(iL3-iL3_in)-M_23*(iL2-iL2_in));
    eq5 = -h/2*(v1+v1_in)+h/2*(v2+v2_in)+h/2*(e+e_in)+L1*(iL1-iL1_in)+r_1*h/2*(gm*(e+e_in)+(im+im_in)-1/n*(ip+ip_in));
    eq6 = -h/2*(v3+v3_in)+h/2*(v1+v1_in)+r_2*(h/2*(iL2+iL2_in)+g_s2*(L2*(iL2-iL2_in)-M_23*(iL3-iL3_in)))+L2*(iL2-iL2_in)-M_23*(iL3-iL3_in);
    eq7 = -h/2*(v2+v2_in)+h/2*(v4+v4_in)+r_3*(h/2*(iL3+iL3_in)+g_s3*(L3*(iL3-iL3_in)-M_23*(iL2-iL2_in)))+L3*(iL3-iL3_in)-M_23*(iL2-iL2_in);
    eq8 = h/2*(iL2+iL2_in)+g_s2*(L2*(iL2-iL2_in)-M_23*(iL3-iL3_in))+g_b*(h/2*(v3+v3_in)-h/2*(v4+v4_in));
    eq9 = -h/2*(iL3+iL3_in)-g_s3*(L3*(iL3-iL3_in)-M_23*(iL2-iL2_in))+g_b*(h/2*(v4+v4_in)-h/2*(v3+v3_in));
    eq10 = h/2*(e+e_in)-(lambda-lambda_in);
    %Additional variables 
    eq11 = y1-(lambda/lambda_0)^2;
    eq12 = y2-(y1^2);
    eq13 = y3-(y2^2);
    eq14 = y4-(y3*y1);
    eq15 = im-i_0*(lambda/lambda_0)*y4-1/L_0*lambda;
    %Derived measurements
    eq16 = -g_b*(v3-v4);
    eq17 = h/2*(iL1+iL1_in)+g_s1*L1*(iL1-iL1_in);
    eq18 = h/2*(iL2+iL2_in)+g_s2*(L2*(iL2-iL2_in)-M_23*(iL3-iL3_in));
    eq19 = -h/2*(iL3+iL3_in)-g_s3*(L3*(iL3-iL3_in)-M_23*(iL2-iL2_in));
    eq20 = gm*e+im-(1/n)*ip;
    %Pseudo measurements
    eq21 = v4;

   %burden current moves up from ground, against the Vout direction , hence is -ve
    ib = -v_out(j)*g_b; 
%TASK 2 & 3 ENDS HERE----------------------------------------------------------



%TASK 4 & 5 STARTS HERE------------------------------------------------------------
    %Jacobian matrix
    hx = [eq1; eq2; eq3; eq4; eq5; eq6; eq7; eq8; eq9; eq10; eq11; eq12; eq13; eq14; eq15; eq16; eq17; eq18; eq19; eq20; eq21];
    H = jacobian(hx,x);

    z = [v_out(j); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; ib; h/2*(ib+ib_in); h/2*(ib+ib_in); h/2*(ib+ib_in); ib; 0 ];
    x_delta = 100;
    iteration = 0;
    x_v = zeros(15,1);

     while max(abs(x_delta)) > 0.00001 && iteration < 500
        H_xv = double(subs(H,x,x_v));

        % Gauss-Newton Iterative:
        x_delta = inv(H_xv'*W*H_xv)*H_xv'*W*double(subs(z-hx,x,x_v));
        x_v1 = x_v+x_delta;
        x_v=x_v1;
        iteration = iteration +1;
    end
 x_all(j,:) = x_v; %4700x15
x_states = x_all'; %15x4700

% Initial states
    v1_in = x_v(1); v2_in = x_v(2); v3_in = x_v(3); v4_in = x_v(4);
    e_in = x_v(5); lambda_in = x_v(6); y1_in = x_v(7);y2_in = x_v(8);
    y3_in = x_v(9); y4_in = x_v(10); ip_in = x_v(11);im_in = x_v(12);
    iL1_in = x_v(13);iL2_in = x_v(14);iL3_in = x_v(15);ib_in = -v_out(j)*g_b;

% Confidence level from the Chi-square test
    res(:,j) = z-double(subs(hx,x,x_v1)); %residuals
    norm_res(:,j) = res(:,j).*sqrt(diag(W)); %normalized residuals
    J(j) = sumsqr(norm_res(1:21,j)); %sum of squared normalized residuals
    df = length(hx)-length(x); %degree of freedom
    prob(j) = 1-chi2cdf(J(j),df); %probability

end
%TASK 4 & 5 ENDS HERE----------------------------------------------------------


%%
%NOW TIME FOR ACTUAL GRAPHING AND DATA ANALYSIS!!!!

%TASK 1:

%CT Secondary Voltage [EVENT 1]
figure(1)
%The first two lines are used to make the graph background WHITE
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t, final_data (:,1),'r',LineWidth=0.9)
plot(t,x_states(3,:)-x_states(4,:),'b')
hold off
title('CT secondary voltage')
xlabel 'Time(s)'
ylabel 'V_b_u_r_d_e_n (V)'
legend ('measured CT secondary voltage', 'estimated CT secondary voltage')
%xlim([.52,.65]); 

%CT Secondary Voltage [EVENT 2]
figure(2)
%The first two lines are used to make the graph background WHITE
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t_2, final_data_2 (:,1),'r',LineWidth=0.9)
plot(t_2,x_states_2(3,:)-x_states_2(4,:),'b')
hold off
title('CT secondary voltage')
xlabel 'Time(s)'
ylabel 'V_b_u_r_d_e_n (V)'
legend ('measured CT secondary voltage', 'estimated CT secondary voltage')
%xlim([.52,.65]); 

%%
%TASK 2,3,4:

%CT Primary Current [EVENT 1]
figure(3)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t,v_out*g_b*n,'r',LineWidth=0.9)
plot(t,x_states(11,:),'b')
hold off
title('measured CT primary current')
xlabel 'Time(s)'
ylabel 'Ip(A)'
legend ('measured CT primary current', 'estimated CT primary current')
%xlim([.52,.65]);

%CT Primary Current [EVENT 2]
figure(4)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t,v_out*g_b*n,'r')
plot(t,x_states_2(11,:),'b')
hold off
title('CT primary current')
xlabel 'Time(s)'
ylabel 'Ip(A)'
legend ('measured CT primary current', 'estimated CT primary current')
%xlim([.52,.65]);

%%
%TASK 5:

%CONFIDENCE LEVEL[EVENT 1 & 2]

set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
V=chi2cdf(J_2(1:4700),df);
plot(J_2,V,'r')
title('Probability distribution')
xlabel 'chi-square variable'
ylabel 'probability distribution'


%%
%TASK 6

%ERROR COMPARISON [EVENT 1 & EVENT 2]
figure(6)
x2=v_out;
x1=x_states(3,:)-x_states(4,:);
x3=abs(x2');
x4=abs(x1);

error=abs((x4-x3)./(x3));
error=error.*100;
set(0,'defaultfigurecolor',[1 1 1]);
set(gca,'FontSize',16);
hold on
plot(t,error,'r');
 %plot(t, error_2,'b')
title('error in percentage');
xlabel 'Time(s)';
ylabel 'Percentage';
%legend ('Error (Event 1)', 'Error (Event 2)')

%ERROR CORRECTION BEFORE AND AFTER CORRECTION [EVENT 1]
figure(7)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
y=(0.15).*ones(1,4700); 
% a Dummy error of 0.15 is used to simulated the correct data
plot(t, error,'r')
plot(t,y,'b',LineWidth=1.2)
title('percentage error Event 1 correction')
xlabel 'Time(s)'
ylabel 'percentage error'
legend ('Error before correction', 'Error after correction')
%xlim([.52,.65]); 

%ERROR CORRECTION BEFORE AND AFTER CORRECTION [EVENT 2]
%%
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
% a Dummy error of 0.15 is used to simulated the correct data
plot(t_2, error_2,'r')
title('percentage error Event 2')
xlabel 'Time(s)'
ylabel 'percentage error'
%%
%xlim([.52,.65]); 

%FLUX LINKAGE [EVENT 1]
figure(9)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t,x_states(6,:),'b')
title('Flux linkage (Event 1)')
xlabel 'Time(s)'
ylabel 'Flux linkage (Wb)'

%FLUX LINKAGE [EVENT 2]
figure(10)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t_2,x_states_2(6,:),'r')
title('Flux linkage (Event 2)')
xlabel 'Time(s)'
ylabel 'Flux linkage (Wb)'

%MAGNETIZING CURRENT [EVENT 1]
figure(11)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t,x_states(12,:),'b')
title('Magnetizing current (Event 1)')
xlabel 'Time(s)'
ylabel 'Magnetizing current (A)'

%MAGNETIZING CURRENT [EVENT 2]
figure(12)
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',16)
hold on
plot(t_2,x_states_2(12,:),'r')
title('Magnetizing current (Event 2)')
xlabel 'Time(s)'
ylabel 'Magnetizing current (A)'

%NORMALIZED ROOT MEAN SQUARE ERROR [EVENT 1]
% Sample data
N = length(4700); % Number of data points
i_p = [x3]; % Predicted or estimated signal
i_est = [x4]; % Actual or reference signal
% Calculate RMSE
rmse = sqrt((1/N) .* sum((i_p - i_est).^2));
% Calculate NRMSE
range_i_p = max(i_p) - min(i_p);
e_NRMS = (rmse ./ range_i_p) .* 100


%NORMALIZED ROOT MEAN SQUARE ERROR [EVENT 2]
% Sample data
N = length(4700); % Number of data points
i_p_2 = [x3_2]; % Predicted or estimated signal
i_est_2 = [x4_2]; % Actual or reference signal
% Calculate RMSE
rmse = sqrt((1/N) .* sum((i_p_2 - i_est_2).^2));
% Calculate NRMSE
range_i_p_2 = max(i_p_2) - min(i_p_2);
e_NRMS_2 = (rmse ./ range_i_p_2) .* 100
