clear; close all; clc

%% Constants and functions
f_c = 100E3;
T = 0.1E-3;
t = 0:1E-7:T;
f = @(t) exp(-1*t.^2);

% Error
E_b = 1;
E_s = @(M) E_b*log2(M);
erfc = @(x) 2/sqrt(pi) * integral(f,x,Inf);
Q = @(x) 1/2*erfc(x/sqrt(2));

% E_b/N_0 between 0 and 10 dB
E_b_N_0_dB = 0:10;
%E_b_N_0 = db2mag(E_b_N_0_dB);
E_b_N_0 = 10.^(E_b_N_0_dB/10);

%% Generate frame of 1024 rand bits
n = 1024;
frame = randi([0 1],1,n);
bit_num = n;

%% 2PSK: Modulate one bit at a time
M = 2;

% Generate each phase shift signal as function handle
s = cell(M,1);
for i = 1:M
    phi = 2*pi*i/M;
    s{i} = @(t) sqrt(2*E_s(M)/T).*cos(2*pi*f_c.*t - phi);
end

% Plot random signal
figure();

for bit = 1:bit_num
    i = frame(bit);

    % Offset time scale by T
    t_offset = t + T*(bit-1);

    % Get signal shift for bit modulation
    signal_func = s{i+1};

    plot(t_offset, signal_func(t_offset));
    hold on;
end
title('2-ary Phase Shfit Key Modulated Signal')
xlabel('t (sec)')
ylabel('signal magnitude')


%% 4PSK: Modulate tow bits at a time
M = 4;
s = cell(M,1);
for i = 1:M
    phi = 2*pi*i/M;
    s{i} = @(t) sqrt(2*E_s(M)/T).*cos(2*pi*f_c.*t-phi);
end

% Plot random signal
figure();
i = 0;
for bit = 1:2:bit_num
    i1 = frame(bit);
    i2 = frame(bit+1);

    t_offset = t + T*i;
    i = i + 1;

    signal_func = s{bi2de([i2 i1]) + 1};

    plot(t_offset, signal_func(t_offset));
    hold on;
end

title(sprintf('%d-ary Phase Shift Key Modulated Signal', M))
xlabel('t (sec)')
ylabel('signal magnitude')

%% 8PSK: Modulate three bits at a time
M = 8;
s = cell(M,1);
for i = 1:M
    phi = 2*pi*i/M;
    s{i} = @(t) sqrt(2*E_s(M)/T).*cos(2*pi*f_c.*t - phi);
end

% Plot random signal
figure();

i = 0;
for bit = 1:3:bit_num - 1
    i1 = frame(bit);
    i2 = frame(bit+1);
    i3 = frame(bit+2);
    t_offset = t + T*i;
    i = i + 1;

    signal_func = s{bit2int([i1 i2 i3]', 3)+1};

    plot(t_offset, signal_func(t_offset));
    hold on;
end

title(sprintf('%d-ary Phase Shift Key Modulated Signal', M));
xlabel('t (sec)')
ylabel('signal magnitude');


%% Errors
% 2 - PSK
P_b2 = zeros(1, length(E_b_N_0));
for dB = 1:11
    x = sqrt(2*E_b_N_0(dB));
    P_b2(dB) = Q(sqrt(2*E_b_N_0(dB)));
end

P_s2 = P_b2;

figure(4);
hold off;
plot(E_b_N_0_dB, P_s2);
hold on;
plot(E_b_N_0_dB, P_b2);
set(gca, 'YScale', 'log')
xlabel('E_{b}/N_{0} (dB)')
ylabel('Error probability')
title('Probabilty of error for 2-PSK')
legend('Symbol error', 'Bit error')

% 4-PSK
P_b4 = P_b2;
P_s4 = zeros(1, length(E_b_N_0));
for dB = 1:11
    P_s4(dB) = erfc(sqrt(E_b_N_0(dB)));
    
end

figure(5);
plot(E_b_N_0_dB, P_s4);
hold on;
plot(E_b_N_0_dB, P_b4);
set(gca, 'YScale', 'log')
xlabel('E_{b}/N_{0} (dB)')
ylabel('Error probability')
title('Probabilty of error for 4-PSK')
legend('Symbol error', 'Bit error')

%% 8-PSK error

P_s8 = zeros(1, length(E_b_N_0));
P_b8 = zeros(1, length(E_b_N_0)); 

for dB = 1:11
    P_s8(dB) = 2 * Q(sqrt(2*E_b_N_0(dB)*log2(8))*sin(pi/8));
    P_b8(dB) = P_s8(dB)/log2(8);
end

figure(6);
plot(E_b_N_0_dB, P_s8);
hold on;
plot(E_b_N_0_dB, P_b8);
set(gca, 'YScale', 'log')
xlabel('E_{b}/N_{0} (dB)')
ylabel('Error probability')
title('Probabilty of error for 8-PSK')
legend('Symbol error', 'Bit error')

%% All error plots
figure(7);

plot(E_b_N_0_dB, P_s2);
hold on;
plot(E_b_N_0_dB, P_b2);
hold on;
plot(E_b_N_0_dB, P_s4);
hold on;
plot(E_b_N_0_dB, P_b4);
hold on;
plot(E_b_N_0_dB, P_s8);
hold on;
plot(E_b_N_0_dB, P_b8);

set(gca, 'YScale', 'log')
xlabel('E_{b}/N_{0} (dB)')
ylabel('Error probability')
title('Probabilty of error for M-PSK')
legend('2-PSK P_{s}', '2-PSK P_{b}','4-PSK P_{s}', '4-PSK P_{b}','8-PSK P_{s}', '8-PSK P_{b}')


