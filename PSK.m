%% Reset Sim
close all; % Close current figures.

%% Define Parameters
M_values = [2, 4, 8];
fc = 100e3;
Eb = 1;
Es = Eb*M_values;
T = 0.1e-3;
fs = 10 * fc;
Ts = 1 / fs;

WGN_SNR = 0;

%% Simulate PSK signals using randomly generated bits
% Number of bits
n = 1024;

% Generate random bits
bits = randi([0, 1], 1, n);

% Anon function to convert bit sequence to M-PSK symbols
% This works by reordering the bits into 2xlog2(M) bit segments and then
% converting it to decimal representations using the reshape and bi2de
% functions respectively.
reshaped_bits = @(bits, M) reshape(bits, log2(M),[]).';
bits_to_symbols = @(bits, M) bi2de(reshaped_bits(bits,M), 'left-msb')';

% Anon function to create M-PSK waveform
psk_waveform = @(symbol, M, fc, Es, T, Ts) sqrt(2 * Es / T) * cos(2 * pi * fc * (0:Ts:T-Ts) - 2 * pi * symbol / M);

% Assign up to 10 colours for plotting
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k', 'r--', 'g--', 'b--'};

% Add loop counter to index Es
figI = 1;
loops = 0;
for M = M_values
    loops = loops + 1;
    % Pad bits with zeros until the length is divisible by log2(M)
    % If we don't do this we get an error
    bits_padded = bits;
    while mod(length(bits_padded), log2(M)) ~= 0
        bits_padded = [bits_padded, 0];
    end
    
    % Convert the bits to M-ary symbols
    symbols = bits_to_symbols(bits_padded, M);

    % Create the M-PSK waveform
    waveform_segments = {};
    noisyWaveformSegments = {};
    for symbol = symbols(1:10)
        waveform_segments{end + 1} = psk_waveform(symbol, M, fc, Es(loops), T, Ts);
        noisyWaveformSegments{end + 1} = awgn(waveform_segments{end}, WGN_SNR);
    end

    % Plot the Normal waveform
    figure(1);
    subplot(1,3,figI)
    hold on;
    time_vector = (0:Ts:(T * 10) - Ts);
    for i = 1:length(waveform_segments)
        start_time = (i - 1) * T;
        plot(start_time + (0:Ts:T-Ts), waveform_segments{i}, colors{symbols(i) + 1});
    end
    hold off;
    title([num2str(M), '-PSK Waveform']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

    % Plot the Noisy waveform
    figure(2);
    subplot(1,3,figI)
    figI = figI + 1;
    hold on;
    time_vector = (0:Ts:(T * 10) - Ts);
    for i = 1:length(noisyWaveformSegments)
        start_time = (i - 1) * T;
        plot(start_time + (0:Ts:T-Ts), noisyWaveformSegments{i}, colors{symbols(i) + 1});
    end
    hold off;
    title([num2str(M), '-PSK Waveform with WGN at a SNR of ',num2str(WGN_SNR),'dB']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

%% Calculate error probabilities
% Define error functions
f = @(t) exp(-1*t.^2);
E_b = 1;
E_s = @(M) E_b*log2(M);
erfc = @(x) 2/sqrt(pi) * integral(f,x,Inf);
Q = @(x) 1/2*erfc(x/sqrt(2));

% E_b/N_0 between 0 and 10 dB
E_b_N_0_dB = 0:10;
E_b_N_0 = 10.^(E_b_N_0_dB/10);

% -- 2-PSK --
P_b2 = zeros(1, length(E_b_N_0));
for dB = 1:11
    x = sqrt(2*E_b_N_0(dB));
    P_b2(dB) = Q(sqrt(2*E_b_N_0(dB)));
end

P_s2 = P_b2;

% -- 4-PSK --
P_b4 = P_b2;
P_s4 = zeros(1, length(E_b_N_0));
for dB = 1:11
    P_s4(dB) = erfc(sqrt(E_b_N_0(dB)));
    
end

% -- 8-PSK --
P_s8 = zeros(1, length(E_b_N_0));
P_b8 = zeros(1, length(E_b_N_0)); 

for dB = 1:11
    P_s8(dB) = 2 * Q(sqrt(2*E_b_N_0(dB)*log2(8))*sin(pi/8));
    P_b8(dB) = P_s8(dB)/log2(8);
end

% -- All error plots --
figure(3);
plot(E_b_N_0_dB, P_s2,'o-.');
hold on;
plot(E_b_N_0_dB, P_b2,'+--');
hold on;
plot(E_b_N_0_dB, P_s4);
hold on;
plot(E_b_N_0_dB, P_b4,'*:');
hold on;
plot(E_b_N_0_dB, P_s8);
hold on;
plot(E_b_N_0_dB, P_b8);

set(gca, 'YScale', 'log')
xlabel('E_{b}/N_{0} (dB)')
ylabel('Error probability')
title('Probabilty of error for M-PSK')
legend('2-PSK P_{s}', '2-PSK P_{b}','4-PSK P_{s}', '4-PSK P_{b}','8-PSK P_{s}', '8-PSK P_{b}')

%% Calculate Spectral efficencies
