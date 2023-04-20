% Parameters
M_values = [2, 4, 8];
fc = 100e3;
Eb = 1;
Es = Eb*M_values;
T = 0.1e-3;
fs = 10 * fc;
Ts = 1 / fs;

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
psk_waveform = @(symbol, M, fc, Eb, T, Ts) sqrt(2 * Eb / T) * cos(2 * pi * fc * (0:Ts:T-Ts) + 2 * pi * symbol / M);

% Assign up to 10 colours for plotting
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k', 'r--', 'g--', 'b--'};

% Add loop counter to index Es
loops = 0;
for M = M_values
    loops = loops + 1;
    % Pad bits with zeros until the length is divisible by log2(M)
    bits_padded = bits;
    while mod(length(bits_padded), log2(M)) ~= 0
        bits_padded = [bits_padded, 0];
    end
    
    % Convert bits to M-ary symbols
    symbols = bits_to_symbols(bits_padded, M);

    % Create M-PSK waveform
    waveform_segments = {};
    for symbol = symbols(1:10)
        waveform_segments{end + 1} = psk_waveform(symbol, M, fc, Es(loops), T, Ts);
    end

    % Plot waveform
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
