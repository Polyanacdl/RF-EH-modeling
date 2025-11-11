clearvars
clc

%% Code implementation for comparison models 
% Authors: [Lacerda, P. C.; Mariano, A. A.; Brante, G.]
% Date: [nov - 11 - 2025]

%% --- BLE Technology - Cmin, Tcharge and Minimal Interval calculation ---

% Parameters extracted from supplementary code from Ref. [9] - https://github.com/imec-idlab/EH-feasibility
Pharv = 1e-6*[80.5:0.5:1e3]; 
%Pharv = 1e-3*[1];
Payload = 20;           
e_sense = 0.05;      
t_sense = 0.01;
v_max = 4.5;
v_min = 2.8;
R = 100000;
 
Ttx_d = (8*(Payload + 18))/R;                       % Ref. [9] (eq.15)
e_sense_tx = e_sense + 6.2e-6 + 16.4e-3*(3*Ttx_d);  % Ref. [9] (eq.16)
Ttx_p = (8*(19 + 18))/R;                            % Ref. [9] (eq.15)
Trd = 0.255;
Trw = 0.255;
e_poll = 28.2e-6 + 16.4e-3*(3*Ttx_p) + 16.1e-6*Trd + 15.6e-3*Trw; % Ref. [9] (eq.17)
pmu_l = 1/1.2;
e_peak = (e_sense_tx + e_poll); % Ref. [9] (eq.18)
i_leak = 10e-6;
p_leak = 45e-6;
t_peak = t_sense + 0.050039 + Ttx_d*3 + Ttx_p*3 + Trd + Trw; % Ref. [9] (eq.19)
cmin = (2*(e_peak/pmu_l + p_leak*t_peak))/(v_max^2 - v_min^2); % Ref. [9] (eq.5)
p_idle = 16.128e-6;
p_cons = p_idle/pmu_l + p_leak; % Ref. [9] (eq.7)
pmu_h = 0.8;
p_harv = pmu_h*Pharv; % Ref. [9] (eq.1)
 
for i = 1:length(Pharv)
t_charge(i) = - ((v_max^2)/p_cons) * cmin * log((v_max - (v_max*p_harv(i))/p_cons)/(v_min - (v_max*p_harv(i))/p_cons)); %equação (8) do paper belga
end
t_charge = t_charge;
 
min_interval = t_charge + t_peak; % Ref. [9] (eq.6)

% Calculations from proposed model's equations
for i = 1:length(Pharv)
teh25(i) = (cmin*(0.25*v_max)^2)/(2*Pharv(i));
taoc25(i) = -teh25(i)/log(1 - (0.25*v_max)/v_max);
teh50(i) = (cmin*(0.50*v_max)^2)/(2*Pharv(i));
taoc50(i) = -teh50(i)/log(1 - (0.50*v_max)/v_max);
teh75(i) = (cmin*(0.75*v_max)^2)/(2*Pharv(i));
taoc75(i) = -teh75(i)/log(1 - (0.75*v_max)/v_max);
taoc_avg(i) = mean([taoc25(i) taoc50(i) taoc75(i)]);
teh_ini(i) = -taoc_avg(i)*log(1 - (0.9999*v_max)/v_max);       % Initial teh 
teh_case1(i) = taoc_avg(i)*log(1 - v_min/v_max) + teh_ini(i);  % Method 1
teh_case2(i) = -taoc_avg(i)*log((v_max - 0.9999*v_max)/(v_max - v_min)); % Method 2 
 
end

charging = teh_case2;

% Note: Resr + R2 = res
res = (v_max - v_min - (i_leak*t_peak)/cmin)/i_leak;   % Equivalent resistance succeeding the supercapacitor
discharging = (v_max - v_min)*cmin/i_leak - cmin*res; % Discharging related to leakage current and equivalent resistance succeeding the supercapacitor
teh = charging + discharging;   % Recharging time due to charging from Vmin to Vmax and discharging due to losses
min_interval_2 = teh + t_peak;  % Minimal interval calculation

% Plot: Comparison between simulated curves from Ref. [9] model and baseline model 
figure(1)
loglog(Pharv*1e6,min_interval,"LineWidth",2)
grid on
hold on
loglog(Pharv*1e6,min_interval_2,"LineWidth",2)
legend('Min. Int. = T_c_h_a_r_g_e + T_p_e_a_k [9]','Min. Int = t_E_H + T_p_e_a_k')
xlabel('Harvested Power [uW]')
ylabel('Minimal Interval [s]')
set(gca,'FontSize',14,'FontWeight','bold')

clearvars teh* taoc* 
clear t_charge min_interval teh_case2 teh_case1 teh teh25 teh50 teh75 taoc25 taoc50 taoc75 taoc_avg teh_ini charging discharging min_interval_2
%% --- LoRaWAN Technology - Cmin, Tcharge and Minimal Interval calculation - SODAQ Explorer off idle ---

% Parameters extracted from supplementary code from Ref. [9] - https://github.com/imec-idlab/EH-feasibility
Pharv = 1e-6*[58:1:1e3];
%Pharv = 1e-3*[0.1];
Payload = 20;
e_sense = 0.05;
t_sense = 0.01;
v_max = 4.5;
v_min = 2.8;
 
% LoRaWAN specification - SODAQ Explorer off idle
pream = 8; 
bw = 125000;
sf = 7;
maxpl = 51; 
cr = 0.8;
ih = 1;
de = 0;
p_tx = 0.134;
e_rx = 0.166;
t_rx = 2.2;
e_setup = 0.456;
t_setup = 13;
p_idle = 0.00000108;
 
s_tx = pream+4.25+8+max(ceil((8*Payload-4*sf+44-20*ih)/(4*(sf-2*de)))*(cr+4),0); % Ref. [9] (eq.9)
t_tx = s_tx*(sf^2/bw); % Ref. [9] (eq.10) 
pmu_l = 1/1.2;
e_peak = (e_sense + e_setup + e_rx + t_tx*p_tx); % Ref. [9] (eq.11)
t_peak = t_setup + t_rx + t_tx + t_sense; % Ref. [9] (eq.12)
i_leak = 10e-6;
p_leak = 45e-6; 
cmin = (2*(e_peak/pmu_l + p_leak*t_peak))/(v_max^2 - v_min^2); % Ref. [9] (eq.5). Note: LoRaWAN without pmu_l
p_cons = p_idle/pmu_l + p_leak; % Ref. [9] (eq.7)
pmu_h = 0.8;
p_harv = pmu_h*Pharv; % Ref. [9] (eq.1)
  
for i = 1:length(p_harv)
t_charge(i) = - ((v_max^2)/p_cons) * cmin * log((v_max - (v_max*p_harv(i))/p_cons)/(v_min - (v_max*p_harv(i))/p_cons)); % Ref. [9] (eq.8)
end
t_charge = t_charge;
 
min_interval = t_charge + t_peak; % Ref. [9] (eq.6)
 
% Calculations from proposed model's equations 
for i=1:length(p_harv)
teh25(i) = (cmin*(0.25*v_max)^2)/(2*Pharv(i));
taoc25(i) = -teh25(i)/log(1 - (0.25*v_max)/v_max);
teh50(i) = (cmin*(0.50*v_max)^2)/(2*Pharv(i));
taoc50(i) = -teh50(i)/log(1 - (0.50*v_max)/v_max);
teh75(i) = (cmin*(0.75*v_max)^2)/(2*Pharv(i));
taoc75(i) = -teh75(i)/log(1 - (0.75*v_max)/v_max);
taoc_avg(i) = mean([taoc25(i) taoc50(i) taoc75(i)]);
teh_ini(i) = -taoc_avg(i)*log(1 - (0.9999*v_max)/v_max);       % Initial teh
teh_case1(i) = taoc_avg(i)*log(1 - v_min/v_max) + teh_ini(i);  % Method 1
teh_case2(i) = -taoc_avg(i)*log((v_max - 0.9999*v_max)/(v_max - v_min)); % Method 2
 
end
 
charging = teh_case2;
 
% Note: Resr + R2 = res
res = (v_max - v_min - (i_leak*t_peak)/cmin)/i_leak;   % Equivalent resistance succeeding the supercapacitor
discharging = (v_max - v_min)*cmin/i_leak - cmin*res; % Discharging related to leakage current and equivalent resistance succeeding the supercapacitor
teh = charging + discharging;   % Recharging time due to charging from Vmin to Vmax and discharging due to losses
min_interval_2 = teh + t_peak;  % Minimal interval calculation

% Plot: Comparison between simulated curves from Ref. [9] model and baseline model
figure(2)
loglog(Pharv*1e6,min_interval,"LineWidth",2)
grid on
hold on
loglog(Pharv*1e6,min_interval_2,"LineWidth",2)
legend('Min. Int. = T_c_h_a_r_g_e + T_p_e_a_k [9]','Min. Int = t_E_H + T_p_e_a_k')
xlabel('Harvested Power [uW]')
ylabel('Minimal Interval [s]')
set(gca,'FontSize',14,'FontWeight','bold')

clearvars teh* taoc* 
clear t_charge min_interval teh_case2 teh_case1 teh teh25 teh50 teh75 taoc25 taoc50 taoc75 taoc_avg teh_ini charging discharging min_interval_2
%% --- 6TiSCH Technology - Cmin, Tcharge and Minimal Interval calculation - CC2538

% Parameters extracted from supplementary code from Ref. [9] - https://github.com/imec-idlab/EH-feasibility
Pharv = 1e-6*[330:1:1e3];
%Pharv = 1e-3*[0.33];
Payload = 20;
e_sense = 0.05;
t_sense = 0.01;
v_max = 4.5;
v_min = 2.8;

nbs = 1;
hops = 1;
t_ts = 0.015;
dao_pl = 96;
eb_pl = 16;
p_tx = 0.000000174;
e_tsoh = 0.0003365;
e_poh = 0.0009706;
p_idle = 0.000180;
overhead = 38;    % hops = 1
e_data = e_sense + e_tsoh + p_tx*(Payload+overhead);
e_dao = e_tsoh + p_tx*p_tx;  
e_peak = (e_poh + max(e_data,e_dao));
if e_data>e_dao
    t_peak = 5*t_ts + t_sense;
else
    t_peak = 5*t_ts + 0;
end

i_leak = 10e-6;
p_leak = 45e-6;
pmu_l = 1/1.2;
cmin = (2*(e_peak/pmu_l + p_leak*t_peak))/(v_max^2 - v_min^2); % Ref. [9] (eq.5)
p_cons = p_idle/pmu_l + p_leak; % Ref. [9] (eq.7)
pmu_h = 0.8;
p_harv = pmu_h*Pharv; % Ref. [7] (eq.1)

for i = 1:length(Pharv)
    t_charge(i) = - ((v_max^2)/p_cons) * cmin * log((v_max - (v_max*p_harv(i))/p_cons)/(v_min - (v_max*p_harv(i))/p_cons)); %equação (8) do paper belga
end
t_charge = t_charge;
 
min_interval = t_charge + t_peak; % Ref. [9] (eq.6)

% Calculations from proposed model's equations
for i = 1:length(Pharv)
teh25(i) = (cmin*(0.25*v_max)^2)/(2*Pharv(i));
taoc25(i) = -teh25(i)/log(1 - (0.25*v_max)/v_max);
teh50(i) = (cmin*(0.50*v_max)^2)/(2*Pharv(i));
taoc50(i) = -teh50(i)/log(1 - (0.50*v_max)/v_max);
teh75(i) = (cmin*(0.75*v_max)^2)/(2*Pharv(i));
taoc75(i) = -teh75(i)/log(1 - (0.75*v_max)/v_max);
taoc_avg(i) = mean([taoc25(i) taoc50(i) taoc75(i)]);
teh_ini(i) = -taoc_avg(i)*log(1 - (0.9999*v_max)/v_max);       % Initial teh
teh_case1(i) = taoc_avg(i)*log(1 - v_min/v_max) + teh_ini(i); % Method 1
teh_case2(i) = -taoc_avg(i)*log((v_max - 0.9999*v_max)/(v_max - v_min)); % Method 2

end

charging = teh_case2;
 
% Note: Resr + R2 = res
res = (v_max - v_min - (i_leak*t_peak)/cmin)/i_leak;   % Equivalent resistance succeeding the supercapacitor
discharging = (v_max - v_min)*cmin/i_leak - cmin*res; % Discharging related to leakage current and equivalent resistance succeeding the supercapacitor
teh = charging + discharging;   % Recharging time due to charging from Vmin to Vmax and discharging due to losses
min_interval_2 = teh + t_peak;  % Minimal interval calculation

% Plot: Comparison between simulated curves from Ref. [9] model and baseline model
figure(3)  
loglog(Pharv*1e6,min_interval,"LineWidth",2) 
grid on 
hold on 
loglog(Pharv*1e6,min_interval_2,"LineWidth",2) 
legend('Min.Int. = T_c_h_a_r_g_e + T_p_e_a_k [9]','Min. Int = t_E_H + T_p_e_a_k') 
xlabel('Harvested Power [uW]') 
ylabel('Minimal Interval [s]') 
set(gca,'FontSize',14,'FontWeight','bold') 

clearvars teh* taoc* 
clear t_charge min_interval teh_case2 teh_case1 teh teh25 teh50 teh75 taoc25 taoc50 taoc75 taoc_avg teh_ini charging discharging min_interval_2
