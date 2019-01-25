function [RE,IM,NP.NB,NT,HDR] = Test_Load_IDEAL

%This loads in phantom IDEAl data for unit testing

DataPath = 'Phantom_data/s_2017012401/';

mtssfpFiles = dir(strcat([DataPath '/' 'mtssfp_' '*.fid']));

for k = mtssfpindex(1):mtssfpindex(2)

    mtssfpFilename = mtssfpFiles(k).name(1:end-4);

    nv = queryPP([experimentPath mtssfpFilename,'.fid'],'nv');

    ne = queryPP([experimentPath mtssfpFilename,'.fid'],'ne');

    te = queryPP([experimentPath mtssfpFilename,'.fid'],'te');

    sw = queryPP([experimentPath mtssfpFilename,'.fid'],'sw');

    pad = queryPP([experimentPath mtssfpFilename,'.fid'],'pad');

    tof = queryPP([experimentPath mtssfpFilename,'.fid'],'tof');

    [RE,IM,NP,NB,NT,HDR] = load_fid([experimentPath mtssfpFilename]);

end


end
