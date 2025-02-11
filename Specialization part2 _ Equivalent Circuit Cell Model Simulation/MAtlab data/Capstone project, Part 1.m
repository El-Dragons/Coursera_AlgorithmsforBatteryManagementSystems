
% Initialize workspace, load the E2 circuit model as well as the E2 dynamic data
addpath readonly
load readonly/E2model.mat; % load parameter values already created for the E2 cell -- this is a single R-C model
% load readonly/E2model2RC.mat; % this is a two R-C model of the E2 cell
load readonly/E2_DYN_P25.mat; % load raw test data for the E2 cell at 25 degC

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (0:deltaT:time(end));
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;

% GRADED FUNCTION (do not modify this line)

% function [vk,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)
%
% ik - current in amperes, where (+) is discharge. Size is N x 1.
% temp  - temperature (degC). Size is N x 1.
% deltaT = sampling interval of data in seconds. Size is 1 x 1 (a scalar)
% model - standard model structure
% z0 - initial SOC. Size is 1 x 1.
% iR0 - initial resistor currents as column vector. Size is Nr x 1 where Nr is 
%       number of R-C pairs in model.
% h0 - initial hysteresis state. Size is 1 x 1.
%
% vest - predicted cell voltage. Size is N x 1.
% rck - predicted resistor currents. Size is N x Nr (first row is set to iR0')
% hk - predicted dynamic hysteresis states. Size is N x 1 (first entry is h0)
% zk - predicted cell state of charge. Size is N x 1 (first entry is z0)
% sik - sign of input current. Size is N x 9
% OCV - predicted cell open circuit voltage. Size is N x 1.
function [vest,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)

  % Force data to be column vector(s) in case user entered data incorrectly
  ik = ik(:); iR0 = iR0(:); temp = temp(:);
  N = length(ik); Nr = length(iR0);
  % initialize some outputs
  vest = zeros(N,1); rck = zeros(N,Nr); hk = zeros(N,1); zk = zeros(N,1); 
  sik = zeros(N,1); OCV = zeros(N,1);
  rck(1,:) = iR0'; hk(1) = h0; zk(1) = z0; sik(1) = 0;
  OCV(1) = OCVfromSOCtemp(z0,temp(1),model);
  
  size(vest)
  size(rck)
  size(hk)
  size(zk)
  size(sik)
  size(OCV)
  % BEGIN MODIFYING CODE AFTER THIS

  T = temp(1:deltaT:length(temp)); % sample code uses only single temperature -- you will need to change this!
  T = T(:);
  % The code reproduced below as a starting point for you is extracted from "simCell.m"
  %
  % Get model parameters from model structure -- notice that these retreive parameter values
  % for only a single temperature. You will need to change this to load parameter values for
  % all temperatures in temp.
  
  Nr
  RCfact(Nr,:) = exp(-deltaT/abs(getParamESC('RCParam',T(1),model)))';
  length(RCfact)
  if length(RCfact) ~= Nr,
    error('iR0 does not have the correct number of entries');
  end

  for k =1:10:length(T)
    G(k) = getParamESC('GParam',T(k),model);
    Q(k) = getParamESC('QParam',T(k),model)
  end
    k = 1
    G(k)= getParamESC('GParam',T(k),model);
    Q(k) = getParamESC('QParam',T(k),model);
    M(k) = getParamESC('MParam',T(k),model);
    M0(k) = getParamESC('M0Param',T(k),model);
    RParam(k) = getParamESC('RParam',T(k),model);
    R0Param(k) = getParamESC('R0Param',T(k),model);
    etaParam(k) = getParamESC('etaParam',T(k),model);
  etaik = ik; etaik(ik<0) = etaParam*ik(ik<0);

  % Simulate the dynamic states of the model
  for k = 2:length(ik),
    rck(k,:) = rck(k-1,:)*diag(RCfact) + (1-RCfact')*etaik(k-1);
    zk(k,:) = zk(k-1,:)-cumsum([0;etaik(k-1,:)])*deltaT/(Q(k)*3600);

  end
  warning('Computing is done');
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  % Hysteresis stuff
  fac=exp(-abs(G*etaik*deltaT/(3600*Q)));
  for k=2:length(ik),
    hk(k)=fac(k-1)*hk(k-1)+(fac(k-1)-1)*sign(ik(k-1));
    sik(k) = sign(ik(k));
    if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
  end
    
  % Compute output equation
  for k=1:length(ik),
    OCV(k,:) = OCVfromSOCtemp(zk(k,:),T(k,:),model);
    vest(k,:) = OCV(k,:) - rck(k,:)*RParam' - R0Param*ik(k,:) + M*hk(k,:) + M0*sik(k,:);
  end 
  
  size(rck)
  size(Q)
  size(zk)
  size(T)
  size(rck)
  size(RParam)
  size(R0Param)
  size(hk)
  size(sik)
  size(M0)
  size(M)
  
  
  
  
  
  % FINISH MODIFYING CODE BEFORE THIS

end
% END GRADED FUNCTION

% Execute simCellTemp to determine voltage and other internal states/variables

%temp = 25*ones(size(current)); % for now, use constant 25 degC temperature.
temp = linspace(25,45,length(current)); % uncomment to simulate temperature ramp 
% Note, you will need to change the default initializations "1,0,0" when using a 2RC model
[vest,rck,hk,zk,sik,OCV] = simCellTemp(current,temp,deltaT,model,1,[1],0);

% The grader will input custom current, temperature, and initial states and compare output
% to expected output. You will not be told what these custom inputs will be, but you can still
% test basic functionality and visualize how the cell behavior changes for different temperatures

% For example, for visualization purposes, plot the measured and simulated voltage data.
% Note that the simulated voltage will not match the measured voltage very well for simulated
% temperatures other than 25 degC since the measured data were collected at 25 degC!

subplot(1,2,1)
plot(time/3600,voltage,time/3600,vest); % factor of 3600 converts seconds -> hours
xlabel('Time (hr)'); ylabel('Voltage (V)'); title('Comparing measured to simulated voltage');
legend('Measured voltage','Simulated voltage');

% Now, plot the voltage prediction error
subplot(1,2,2)
plot(time/3600,1000*(voltage(:)-vest(:)));
xlabel('Time (hr)'); ylabel('Voltage (mV)'); title('Voltage prediction error');
