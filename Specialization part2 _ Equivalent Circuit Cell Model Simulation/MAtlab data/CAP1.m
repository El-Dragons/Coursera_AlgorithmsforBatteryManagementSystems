% Initialize workspace, load the E2 circuit model as well as the E2 dynamic data
% addpath readonly

% Resample at consistent 1Hz rate.
deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (0:deltaT:time(end));
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;
temp = 25*ones(size(current)); 
temp = linspace(25,45,length(current));
temp = linspace(25,45,length(current));
  ik = current; 
  z0 =1;
  h0=0.2;
  iR0 = 0.1;

  % Force data to be column vector(s) in case user entered data incorrectly
  ik = ik(:); 
  iR0 = iR0(:); 
  temp = temp(:);
  N = length(ik); 
  Nr = length(iR0);
  % initialize some outputs
  vest = zeros(N,1); 
  rck  = zeros(N,Nr); 
  hk   = zeros(N,1); 
  zk   = zeros(N,1); 
  sik  = zeros(N,1); 
  OCV  = zeros(N,1);
  
  rck(1,:) = iR0'; 
  hk(1)    = h0; 
  zk(1)    = z0; 
  sik(1)   = 0;
  OCV(1)   = OCVfromSOCtemp(z0,temp(1),model);

  % BEGIN MODIFYING CODE AFTER THIS
  T = temp; % sample code uses only single temperature -- you will need to change this!
  % The code reproduced below as a starting point for you is extracted from "simCell.m"
  %
  % Get model parameters from model structure -- notice that these retreive parameter values
  % for only a single temperature. You will need to change this to load parameter values for
  % all temperatures in temp.
  %warning(" checkpoint 0 ! ")
  
  RCfact(1,:) = exp(-deltaT./abs(getParamESC('RCParam',temp(1),model)))';
  
  %warning(" checkpoint 0 Ended ! ")
  
  if length(RCfact) ~= Nr,
    error('iR0 does not have the correct number of entries');
  end
  
  %warning(" checkpoint 1 ! ")
  etaik    = ik(:);
  G        = zeros(N,1);
  Q        = zeros(N,1);
  etaParam = zeros(N,1);
  M        = zeros(N,1);
  M0       = zeros(N,1);
  R0Paramm = zeros(N,1);
  RParam   = zeros(N,Nr);      
  for i=1:1:N,
    M(i)       = getParamESC('MParam',T(i),model);
    M0(i)      = getParamESC('M0Param',T(i),model);
    R0Param(i) = getParamESC('R0Param',T(i),model);
    RCfact(i,:)= exp(-deltaT./abs(getParamESC('RCParam',T(i),model)))';
    RParam(i,:)= getParamESC('RParam',T(i),model);
    etaParam(i)= getParamESC('etaParam',T(i),model);
    Q(i)       = getParamESC('QParam',T(i),model);
    G(i)       = getParamESC('GParam',T(i),model);
    
    if ik(i) < 0,
      etaik(i) = etaParam(i)*ik(i);
    end
  
  end
