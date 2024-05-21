clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% 讀 SorPos.mat 檔 (SorPos) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorPosNum = 4;                                           % number of source positions
Sor_spacing = 0.3;                                       % source spacing
reverberation_time = 0.6;                                % Reverberation time (s)
points_rir = 12288;                                      % Number of rir points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SorPos_filename_str = ['h\SorPos_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), '.mat'];
SorPos_filemane = join(SorPos_filename_str, '');
load(SorPos_filemane)

%% 讀 h.mat 檔  (h) %%
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

%% 讀音檔 (source) %%
fs = 16000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*fs;

Second_total = Second*SorPosNum;
SorLen_total = SorLen*SorPosNum;

[source_transpose, fs] = audioread('245.wav', [1, SorLen_total]);    % speech source
source = source_transpose.';

%% RIR mix source (y_delay and y_nodelay) %%
% convolution source and RIR %
MicNum = 30;
as = zeros(MicNum, points_rir+SorLen-1, SorPosNum);
for j = 1:SorPosNum 
    for i = 1:MicNum
        as(i, :, j) = conv(h(i, :, j), source(:, 1+(j-1)*SorLen:j*SorLen));
    end
end

% introduce delay %
NFFT = 1024;
hopsize = 256;
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF

% generate y_delay and y_nodelay %
y_delay = zeros(MicNum, SorLen_total);
y_nodelay = zeros(MicNum, SorLen_total);
for j = 1:SorPosNum
    y_delay(:, 1+extra_delay_y+(j-1)*SorLen:j*SorLen) = as(:, 1:SorLen-extra_delay_y, j);
    y_nodelay(:, 1+(j-1)*SorLen:j*SorLen) = as(:, 1:SorLen, j);
end

% 存 y_delay and y_nodelay .mat 檔 %
y_delay_filename_str = ['y\y_delay_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.mat'];
y_delay_filename = join(y_delay_filename_str, '');
save(y_delay_filename, 'y_delay')

y_nodelay_filename_str = ['y\y_nodelay_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.mat'];
y_nodelay_filename = join(y_nodelay_filename_str, '');
save(y_nodelay_filename, 'y_nodelay')

%% WPE (y_wpe) %%
% do wpe %
y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
y_wpe = y_wpe.';

% 存 wpe .mat 檔 %
y_wpe_filename_str = ['y\y_wpe_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), 'x', string(Second), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
save(y_wpe_filename, 'y_wpe')
