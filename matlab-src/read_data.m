% addpath
addpath('../dataset/NTU_wind/');
filename = 'marine_dr_complete.csv';

data = importdata(filename);
time = datestr(data.textdata);
ts1 = timeseries(data.data(:,1), time, 'name', 'speed');%data.textdata{1,2});
ts1.datainfo.units='m/s';
ts2 = timeseries(data.data(:,2), time, 'name', data.textdata{1,3});
ts2.datainfo.units='degree';
% ts3 = timeseries(data.data(:,3), time, 'name', data.textdata{1,4});
% ts3.datainfo.units='m/s';
% ts4 = timeseries(data.data(:,4), time, 'name', data.textdata{1,5});
% ts4.datainfo.units='degree';
tsc = tscollection({ts1 ts2}, 'name', 'Marine Dr Wind');
tsc.timeinfo.units='10mins';
tsc.timeinfo.format='ddmmmyyyy hh:mm:ss';

ts1_summary=[min(ts1.data), quantile(ts1.data,0.25), median(ts1.data), quantile(ts1.data, 0.75), max(ts1.data)];

