function HRV = ecg_extraction(data,fs)
%% This function is used to extract features of ECG
% Input:
%   data:   ECG signal
%   fs:     sample frequency
%
% Output:
%   R_value, R_loc: value and location of R peak
%   Q_value, Q_loc: value and location of Q peak
%   S_value, S_loc: value and location of S peak
%   J_value, J_loc: value and location of J peak
%   K_value, K_loc: value and location of K peak
%   T_value, T_loc: value and location of T peak
%   P_value, P_loc: value and location of P peak
%   RR:             RR interval
%   PR:             PR segment
%   QT:             QT segment
%   BPM:            Beat per min
%   tqrs:           duration of QRS complex
%   trr:            duration of RR interval
%   tpr:            duration of PR segment
%   tqt:            duration of QT segment

d = data;
Fs = fs;
% t = [0:length(d)-1]/fs;

if isrow(d) == 0
    d = d';
end

%% ============================ filter baseline ===========================
[a,b] = butter(5,[0.5 30]/(Fs/2));
c = filtfilt(a,b,d);
c = c / max(c);

%% ============================ low and high pass filter =========================== 
%LPF
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a = [1 -2 1];
h_LP = filter(b,a,[1 zeros(1,12)]);

x2 = conv (c ,h_LP);
x2 = x2 (6+(1: length(d))); %cancle delay
x2 = x2 / max(x2);

%HPF
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];

h_HP = filter(b,a,[1 zeros(1,32)]); % impulse response of HPF
 
x3 = conv (x2 ,h_HP);
x3 = x3 (16+(1: length(d))); %cancle delay
x3 = x3 / max(x3);

%% ======================== Make impulse response =========================
h = [-1 -2 0 2 1]/8;

% Apply filter
x4 = conv (x3 ,h);
x4 = x4 (2+(1: length(d)));
x4 = x4 / max(x4);

%% ============================== Squaring ================================ 
x5 = x4 .^2;
x5 = x5/ max(x5);

%% ====================== Moving-window integrator ======================== 
n = 31;
h = ones (1 ,n)/n;
% Delay = 15; % Delay in samples

% Apply filter
x6 = conv (x5 ,h);
x6 = x6 (15+(1: length(d)));
x6 = x6 / max(x6);

max_h = max(x6);
thresh = mean (x6);
poss_reg =(x6>thresh*max_h);

left = find(diff([0 poss_reg])==1);
right = find(diff([poss_reg 0])==-1);


%% ============= R peak detection and rejection the fail peak =============
% ====== detect R peak ======
R_value = zeros(1,length(left));
R_loc = zeros(1,length(left));

for i=1:length(left)
    [R_value(i), R_loc(i)] = max(c(left(i):right(i)) );
    R_loc(i) = R_loc(i) - 1 + left(i); % add offset
end

need_to_remove = zeros(1,length(left));
j = 1;
for i = 3:length(R_loc) - 2
    value1 = (R_loc(i) - R_loc(i - 1)) / (R_loc(i - 1) - R_loc(i - 2));
    value2 = (R_loc(i + 1) - R_loc(i - 1)) / (R_loc(i - 1) - R_loc(i - 2));
    value3 = (R_loc(i + 2) - R_loc(i + 1)) / (R_loc(i - 1) - R_loc(i - 2));
    if value1 < 0.5 && abs(value2 - value3) < 0.75
        need_to_remove(j) = i;
        j = j + 1;
    end;
end;
need_to_remove(need_to_remove == 0) = [];

R_loc(need_to_remove) = [];
% R_value(need_to_remove) = [];
HRV = one(1,length(R_loc)-1);
RR = zeros(1,length(R_loc)-1);
trr = zeros(1,length(R_loc)-1);

%% =================== detect other peaks and interval ====================
for i = 2:length(R_loc)   
    % ====== RR interval ====== 
    RR(i-1) = R_loc(i)-R_loc(i-1);
    trr(i-1) = RR(i-1)/Fs;

    % ====== BPM (vent rate) ====== 
    HRV(i-1) = 60/trr(i-1);
end


