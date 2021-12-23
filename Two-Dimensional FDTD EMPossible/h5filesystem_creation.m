% INITIALIZE MATLAB
close all;
clc;
clear all;

% TIME STAMP
t1 = clock;

% Import dataset

load_time1 = clock;
disp('Loading the dataset...');
load('sample_set_50.mat');
load_time2 = clock;
load_time = etime(load_time2, load_time1);
disp(['Elapsed time is ' num2str(load_time) ' seconds.']);
disp(['Elapsed time is ' num2str(load_time/60) ' minutes.']);
disp('Dataset Loaded!');

% Create hdf5 file
filename = 'sample_comparison_50.h5';

[rows_mat, cols_mat] = size(Spectrum(1).DevER);
[rows_opt, cols_opt] = size(Spectrum(1).A);

for index = 1:length(Spectrum)
    h5create(filename,strcat('/material/DevER/', num2str(index-1)),[rows_mat cols_mat], 'ChunkSize',[20 20],'Deflate',9)
    h5create(filename,strcat('/material/DevSIG/', num2str(index-1)),[rows_mat cols_mat], 'ChunkSize',[20 20],'Deflate',9)
    h5create(filename,strcat('/optical/A/', num2str(index-1)) ,[rows_opt cols_opt], 'ChunkSize',[1 20],'Deflate',9)
    h5create(filename,strcat('/optical/R/', num2str(index-1)) ,[rows_opt cols_opt], 'ChunkSize',[1 20],'Deflate',9)
    h5create(filename,strcat('/optical/T/', num2str(index-1)) ,[rows_opt cols_opt], 'ChunkSize',[1 20],'Deflate',9)
    h5create(filename,strcat('/optical/C/', num2str(index-1)) ,[rows_opt cols_opt], 'ChunkSize',[1 20],'Deflate',9)
end

for file_index = 1:length(Spectrum)
    h5write(filename,strcat('/material/DevER/', num2str(file_index-1)), Spectrum(file_index).DevER);
    h5write(filename,strcat('/material/DevSIG/', num2str(file_index-1)), Spectrum(file_index).DevSIG);
    h5write(filename,strcat('/optical/A/', num2str(file_index-1)), Spectrum(file_index).A);
    h5write(filename,strcat('/optical/R/', num2str(file_index-1)), Spectrum(file_index).R);
    h5write(filename,strcat('/optical/T/', num2str(file_index-1)), Spectrum(file_index).T);
    h5write(filename,strcat('/optical/C/', num2str(file_index-1)), Spectrum(file_index).C);
end
h5disp(filename);


t2 = clock;
t = etime(t2,t1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);