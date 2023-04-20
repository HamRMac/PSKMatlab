%% Define Parameters
M_values = [2, 4, 8];
fc = 100e3;
Eb = 1;
Es = Eb*M_values;
T = 0.1e-3;
fs = 10 * fc;
Ts = 1 / fs;

%% Simulate PSK signals using randomly generated bits
% Number of bits
n = 1024;

% Generate random bits
bits = randi([0, 1], 1, n);

% Anon function to convert bit sequence to M-PSK symbols
% This works by reordering the bits into 2xlog2(M) bit segments and then
% converting it to decimal representations using the reshape and bi2de
% functions respectively.
reshaped_bits = @(bits, M) reshape(bits, [], log2(M));
bits_to_symbols = @(bits, M) bi2de(reshaped_bits(bits,M), 'left-msb')';

% Anon function to create M-PSK waveform
psk_waveform = @(symbol, M, fc, Es, T, Ts) sqrt(2 * Es / T) * cos(2 * pi * fc * (0:Ts:T-Ts) - 2 * pi * symbol / M);

% Assign up to 10 colours for plotting
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k', 'r--', 'g--', 'b--'};

% Add loop counter to index Es
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
    for symbol = symbols(1:10)
        waveform_segments{end + 1} = psk_waveform(symbol, M, fc, Es(loops), T, Ts);
    end

    % Plot the waveform
    figure;
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
end

%% Calculate error probabilities
% Define error functions
syms erfc_probFunc(t) erfc(x) Q(x)
erfc_probFunc(t) = exp(-t.^2);
erfc(x) = 2/sqrt(pi)*int(erfc_probFunc,x,Inf);
Q(x) = (1/2)*erfc(x/sqrt(2));

% Define Parameters
Eb_N0 = 10.^((0:10)/10);

% 2PSK error
error_2PSK_bit = Q(sqrt(2.*Eb_N0));
error_2PSK_sym = error_2PSK_bit;

% 4PSK error
error_4PSK_bit = Q(sqrt(2.*Eb_N0));
error_4PSK_sym = erfc(sqrt(Eb_N0));

% 8PSK error
error_8PSK_sym = 2.*Q(sqrt(2.*Eb_N0.*log2(8)).*sin(pi/8));
error_8PSK_bit = error_8PSK_sym/log2(8);
% plotting somethign
figure
hold on
loglog(Eb_N0,error_2PSK_bit,'r-')
loglog(Eb_N0,error_2PSK_sym,'r--')
loglog(Eb_N0,error_4PSK_bit,'g-')
loglog(Eb_N0,error_4PSK_sym,'g--')
loglog(Eb_N0,error_8PSK_bit,'b-')
loglog(Eb_N0,error_8PSK_sym,'b--')
title('Plot of M-PSK error probabilities for bits and symbols');
xlabel('Eb/N0 (SNR)');
ylabel('Probability of errors');
grid on;