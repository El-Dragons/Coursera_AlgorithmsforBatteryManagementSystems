load E2model.mat; % load parameter values already created for the E2 cell -- this is a single R-C model
 %load readonly/E2model2RC.mat; % this is a two R-C model of the E2 cell
load E2_DYN_P25.mat; % load raw test data for the E2 cell at 25 degC

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (0:deltaT:time(end));
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;