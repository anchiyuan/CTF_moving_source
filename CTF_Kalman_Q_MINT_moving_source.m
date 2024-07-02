clc; clear;
close all;

%% load data and set basic parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorPosNum = 3;                                           % number of source positions
Sor_spacing = 0.1;                                       % source spacing
reverberation_time = 0.4;                                % Reverberation time (s)
points_rir = 8192;                                      % Number of rir points
look_mic = 38;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load SorPos %
SorPos_filename_str = ['h\SorPos_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), '.mat'];
SorPos_filemane = join(SorPos_filename_str, '');
load(SorPos_filemane)

% load h %
for j = 1:SorPosNum
    rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), '.', string(j), '.mat'];
    rir_filemane = join(rir_filename_str, '');
    if j == 1
        h = cell2mat(struct2cell(load(rir_filemane)));
    else
        h_load = cell2mat(struct2cell(load(rir_filemane)));
        h(:, :, j) = h_load;
    end

end

% load source %
fs = 16000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*fs;

Second_total = Second*SorPosNum;
SorLen_total = SorLen*SorPosNum;

[source_transpose, fs] = audioread('245.wav', [1, SorLen_total]);    % speech source
source = source_transpose.';

% load y_delay and y_nodelay %
y_delay_filename_str = ['y\y_delay_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.mat'];
y_delay_filename = join(y_delay_filename_str, '');
load(y_delay_filename)

y_nodelay_filename_str = ['y\y_nodelay_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.mat'];
y_nodelay_filename = join(y_nodelay_filename_str, '');
load(y_nodelay_filename)

% load y_wpe %
y_wpe_filename_str = ['y\y_wpe_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

%% STFT of all needed data %%
% windows %
NFFT = 1024;
hopsize = 256;
win = hamming(NFFT);
osfac = round(NFFT/hopsize);
frequency = NFFT/2 + 1;
freqs_vector = linspace(0, fs/2, frequency);
L = length(hopsize:hopsize:points_rir+2*NFFT-2);    % (len(win) + len(win) - 1) + points_rir - 1

% STFT of y_delay (Y_delay) %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

% STFT of y_wpe (Y_wpe) %
y_wpe_transpose = y_wpe.';
[Y_wpe, ~, ~] = stft(y_wpe_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% specify frame index %%
NumOfFrame = size(Y_wpe, 2);

frame_index = zeros(SorPosNum, 2);
for j = 1:SorPosNum
    if j == 1
        frame_index(1, 1) = 1;
        frame_index(1, 2) = ceil((SorLen*j-NFFT)/hopsize) + 1;
    else
        frame_index(j, 1) = frame_index(j-1, 2) + 1;
        frame_index(j, 2) = ceil((SorLen*j-NFFT)/hopsize) + 1;
    end

end

frame_index(SorPosNum, 2) = NumOfFrame;

%% DAS beamformer %%
MicNum = 38;
SorNum = 1;    % 與 SorPosNum 沒關係 為了 Y_DAS 而設的
c = 343;

% distributed 8 mic %
mic_x = [ 200 ; 300 ; 300 ; 200 ; 200 ; 300 ; 300 ; 200 ]./100;
mic_y = [ 200 ; 200 ; 300 ; 300 ; 200 ; 200 ; 300 ; 300 ]./100;
mic_z = [ 100 ; 100 ; 100 ; 100 ; 200 ; 200 ; 200 ; 200 ]./100;
MicPos = [mic_x, mic_y, mic_z,];

% ULA 30 mics %
MicNum_TDOA = 8;  
MicStart = [210, 200, 100]/100;
spacing = 0.02;
for i = MicNum_TDOA+1:MicNum
    MicPos(i, :) = [MicStart(1, 1)+(i-(MicNum_TDOA+1))*spacing, MicStart(1, 2), MicStart(1, 3)];
end

Y_DAS = zeros(frequency, NumOfFrame);
for j = 1:SorPosNum
    % 算 mic 與 source 之距離 %
    distance = zeros(MicNum, SorNum);
    for i = 1 : MicNum
        distance(i, :) =  sqrt(sum((SorPos(j, :) - MicPos(i, :)).^2));
    end
    
    % 算 a %
    a = zeros(MicNum, SorNum, frequency);
    for n = 1:frequency
        omega = 2*pi*freqs_vector(n);
        a(:, :, n) = exp(-1j*omega/c*distance)./distance;
    end
    
    % 算 DAS weight %
    w = a/MicNum;
    
    % 算 Y_DAS %
    for FrameNo= frame_index(j, 1):frame_index(j, 2)
        for n = 1: frequency
             Y_DAS(n, FrameNo) = w(:, :, n)'*squeeze(Y_wpe(n, FrameNo, :));
        end  
    
    end

end

%% predict CTF with Kalman nonstationary filter (A) %%
start_frame = L;
end_frame = NumOfFrame;

% Kalman nonstationary filter %
A = zeros(MicNum, L, frequency, SorPosNum);
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        Q = 10^(-4)*eye(L);    % process noise covariance matrix %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_frame:end_frame
            % time update %
            P = P + Q;

            % measurement update %
            K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;

            % 轉換的前一個 frame 紀錄 weight %
            for j = 1:SorPosNum
                if FrameNo == frame_index(j, 2)-1
                    A(i, :, n, j) = weight';
                end

            end

        end
    
    end

end

toc

%% A 轉回時域 (A_tdomain) %%
NRMSPM = zeros(SorPosNum, 1);

f1 = figure;
f2 = figure;

for j = 1:SorPosNum
    A_forplot = zeros(frequency, L, MicNum);
    for i = 1 : MicNum
        A_forplot(:, :, i) = squeeze(A(i, :, :, j)).';
    end
    
    A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
    A_tdomain = A_tdomain(:, hopsize*(osfac-1)+1:end);
    ratio_A_tdomain = zeros(MicNum, 1);
    for i = 1:MicNum
        ratio_A_tdomain(i, :) = max(abs(h(i, :, j)))/max(abs(A_tdomain(i, :)));
    end
    
    A_tdomain = A_tdomain.*ratio_A_tdomain;

    ATF = fft(h(:, :, j), points_rir, 2);
    ATF_estimated = fft(A_tdomain, points_rir, 2);

    % 畫 A_tdomain time plot %
    figure(f1)
    subplot(SorPosNum, 1, j);
    plot(h(look_mic, :, j), 'r');
    hold on
    plot(A_tdomain(look_mic, :), 'b');
    hold off
    h_yaxis_upperlimit = max(h(look_mic, :, j)) + 0.01;
    h_yaxis_underlimit = min(h(look_mic, :, j)) - 0.01;
    ylim([h_yaxis_underlimit h_yaxis_upperlimit])
    title(['position = ',num2str(j)])
    legend('ground-truth RIR', 'estimated RIR')
    xlabel('points')
    ylabel('amplitude')

    figure(f2)
    subplot(SorPosNum, 2, j*2-1);
    semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF(look_mic, 1:points_rir/2+1))), 'r');
    hold on
    semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b');
    hold off
    xlim([200 8000])
    legend('ground-truth ATF', 'estimated ATF')
    xlabel('frequency (Hz)')
    ylabel('dB')
    title(['magnitude  position = ',num2str(j)])

    subplot(SorPosNum, 2, j*2);
    semilogx(linspace(0, fs/2, points_rir/2+1), unwrap(angle(ATF(look_mic, 1:points_rir/2+1))), 'r');
    hold on
    semilogx(linspace(0, fs/2, points_rir/2+1), unwrap(angle(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b');
    hold off
    xlim([200 8000])
    legend('ground-truth ATF', 'estimated ATF')
    xlabel('frequency (Hz)')
    ylabel('phase (radius)')
    title(['phase  position = ',num2str(j)])

    % 算 NRMSPM %
    h_NRMSPM = reshape(h(:, :, j).', [MicNum*points_rir 1]);
    aa_NRMSPM = reshape(A_tdomain.', [MicNum*points_rir 1]);
    NRMSPM(j, :) = 20*log10(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));
end

% save fig %
fig_filename_str = ['fig\RIR_nonstationary_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.fig'];
fig_filename = join(fig_filename_str, '');
saveas(f1, fig_filename)

fig_filename_str = ['fig\ATF_nonstationary_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.fig'];
fig_filename = join(fig_filename_str, '');
saveas(f2, fig_filename)

fprintf('done\n')
