clear all
clc
%% Get DGPs
dgps = dir("Data");
dgps = dgps([dgps.isdir]);
dgps = dgps(~ismember({dgps.name}, {'.', '..'}));
%% Get files in DGP
for i = 1:length(dgps)
    if ~exist(fullfile("Results",dgps(i).name),'dir')
        mkdir(fullfile("Results",dgps(i).name));
    end
    files = dir(fullfile("Data",dgps(i).name));
    files = files(~ismember({files.name}, {'.', '..'}));
    for j = 1:length(files)
        sample = fullfile("Data",dgps(i).name,files(j).name);
        res_file = strrep(sample,"mat","csv");
        res_file = strrep(res_file,"Data","Results");
        run_sims;
    end
end