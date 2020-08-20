%Final Capston 1 Project.m
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
% vest - predicted cell voltage. Size is N x 1.
% rck - predicted resistor currents. Size is N x Nr (first row is set to iR0').
% hk - predicted dynamic hysteresis states. Size is N x 1 (first entry is h0).
% zk - predicted cell state of charge. Size is N x 1 (first entry is z0).
% sik - sign of input current. Size is N x 1.
% OCV - predicted cell open circuit voltage. Size is N x 1.

function [vest,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)

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
  T = temp(:); % sample code uses only single temperature -- you will need to change this!
  % The code reproduced below as a starting point for you is extracted from "simCell.m"
  %
  % Get model parameters from model structure -- notice that these retreive parameter values
  % for only a single temperature. You will need to change this to load parameter values for
  % all temperatures in temp.
  %warning(" checkpoint 0 ! ")
  
  RCfact(1,:) = exp(-deltaT./abs(getParamESC('RCParam',T(1),model)))';
  
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
   
  %warning("Checkpoint 2 !")
  % etaik(ik<0) = etaParam(ik<0)*ik(ik<0);
  % Simulate the dynamic states of the model
  for k = 2:length(ik),
    present = rck(k-1,:)*diag(RCfact(k-1,:)) +(1-RCfact(k-1,:))*etaik(k-1);
    rck(k,:) = present;
  end

  %warning("Checkpoint 3 !")
  for k = 2:length(ik),
    zk(k) = zk(k-1) - etaik(k-1)*deltaT/(Q(k)*3600);
  end
 % cumsum(1) = 0;
 % for k = 2:length(ik), 
 %   zk(k) = zk(1) -cumsum(k-1);% etaik(k)*deltaT/(Q(k)*3600);
 %   cumsum(k) =cumsum(k-1) + etaik(k-1)*deltaT/(Q(k-1)*3600);
 % end
  %zk = z0-cumsum([0;etaik(1:end-1)])*deltaT/(Q*3600);
  
  %warning("checkpoint 4 !")
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  % Hysteresis stuff
  for i= 1:1:length(ik),
   fac(i)=exp(-abs(G(i)*etaik(i)*deltaT/(3600*Q(i)))); 
  end
  %fac=exp(-abs(G*etaik.*deltaT/(3600*Q)));

  %warning("checkpoint 5 !!")
  for k=2:length(ik),
    hk(k)  = fac(k-1)*hk(k-1)+(fac(k-1)-1)*sign(ik(k-1));
    sik(k) = sign(ik(k));
    if abs(ik(k))<Q(k)/100, 
      sik(k) = sik(k-1); 
    end
  end
   %warning("checkpoint 6 !!")  
  % Compute output equation
  for k=1:1:N,
    OCV(k) = OCVfromSOCtemp(zk(k),T(k),model);
  end
  
  % warning("checkpoint 7 !!")
  %vest = OCV - rck(1,:).*RParam(1,:) - R0Param.*ik + M.*hk + M0.*sik;
  for k=1:1:N,
    vest(k) = OCV(k) - rck(k,Nr)*RParam(k,Nr)- R0Param(k)*ik(k) + M(k)*hk(k) + M0(k)*sik(k);
  end
  %warning('End of Execution !!')
  %vest = OCV - rck.*RParam' - R0Param.*ik + M.*hk + M0.*sik;
  % FINISH MODIFYING CODE BEFORE THIS
end
% END GRADED FUNCTION