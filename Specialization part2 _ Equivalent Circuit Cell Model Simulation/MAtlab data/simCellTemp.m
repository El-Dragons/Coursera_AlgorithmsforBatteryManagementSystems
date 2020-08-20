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

  % BEGIN MODIFYING CODE AFTER THIS

  T =temp(:); % sample code uses only single temperature -- you will need to change this!
  
  % The code reproduced below as a starting point for you is extracted from "simCell.m"
  %
  % Get model parameters from model structure -- notice that these retreive parameter values
  % for only a single temperature. You will need to change this to load parameter values for
  % all temperatures in temp.
  
  for k =1:Nr
      RCfact(k,1)= exp(-deltaT./abs(getParamESC('RCParam',temp(1,:),model)))';
  end
  if length(RCfact) ~= Nr,
    error('iR0 does not have the correct number of entries');
  end
  G = zeros(size(T));Q = zeros(size(T));M = zeros(size(T));M = zeros(size(T)); M0 = zeros(size(T));
  RParam = zeros(size(T));R0Param = zeros(size(T));etaParam = zeros(size(T));
  for k=1:N,
%   k=1;
    G(k) = getParamESC('GParam',T(k),model); 
    Q(k) = getParamESC('QParam',T(k),model);
    etaParam(k) = getParamESC('etaParam',T(k),model);
    M(k) = getParamESC('MParam',T(k),model);
    M0(k) = getParamESC('M0Param',T(k),model);
    RParam(k) = getParamESC('RParam',T(k),model);
    R0Param(k) = getParamESC('R0Param',T(k),model);
  end
  etaik = ik; etaik(ik<0) = etaParam(ik<0).*ik(ik<0);

  % Simulate the dynamic states of the model
  for k = 2:length(ik),
    rck(k,:) = rck(k-1,:)*diag(RCfact) + (1-RCfact')*etaik(k-1);
     %rck(k,:) = rck(k-1,:)*RCfact(:,1)+ (1-RCfact(N:,1))*etaik(k-1);
  end

  for k =2:N,
      zk(k) = zk(k-1)-etaik(k)*deltaT/(Q(k)*3600);
  end
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  % Hysteresis stuff
  for k=1:N
    fac(k)=exp(-abs(G(k)*etaik(k)*deltaT/(3600*Q(k))));
  end
  
  for k=2:length(ik),
    hk(k)=fac(k-1)*hk(k-1)+(fac(k-1)-1)*sign(ik(k-1));
    sik(k,:) = sign(ik(k));
    if abs(ik(k))<Q(k)/100, sik(k) = sik(k-1); end
  end
    
  % Compute output equation
  for k= 1:N
    OCV(k) = OCVfromSOCtemp(zk(k),T(1),model);
  end
  for k=1:N
    vest(k,1) = OCV(k,1) - rck(k,1)*RParam(k) - R0Param(k)*ik(k,1) + M(k)*hk(k,1) + M0(k)*sik(k,:);
  end
  % FINISH MODIFYING CODE BEFORE THIS

end









% % GRADED FUNCTION (do not modify this line)
% 
% % function [vk,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)
% %
% % ik - current in amperes, where (+) is discharge. Size is N x 1.
% % temp  - temperature (degC). Size is N x 1.
% % deltaT = sampling interval of data in seconds. Size is 1 x 1 (a scalar)
% % model - standard model structure
% % z0 - initial SOC. Size is 1 x 1.
% % iR0 - initial resistor currents as column vector. Size is Nr x 1 where Nr is 
% %       number of R-C pairs in model.
% % h0 - initial hysteresis state. Size is 1 x 1.
% %
% % vest - predicted cell voltage. Size is N x 1.
% % rck - predicted resistor currents. Size is N x Nr (first row is set to iR0')
% % hk - predicted dynamic hysteresis states. Size is N x 1 (first entry is h0)
% % zk - predicted cell state of charge. Size is N x 1 (first entry is z0)
% % sik - sign of input current. Size is N x 9
% % OCV - predicted cell open circuit voltage. Size is N x 1.
% function [vest,rck,hk,zk,sik,OCV] = simCellTemp(ik,temp,deltaT,model,z0,iR0,h0)
% 
%   % Force data to be column vector(s) in case user entered data incorrectly
%   ik = ik(:); iR0 = iR0(:); temp = temp(:);
%   N = length(ik); Nr = length(iR0);
%   % initialize some outputs
%   vest = zeros(N,1); rck = zeros(N,Nr); hk = zeros(N,1); zk = zeros(N,1); 
%   sik = zeros(N,1); OCV = zeros(N,1);
%   rck(1,:) = iR0'; hk(1) = h0; zk(1) = z0; sik(1) = 0;
%   OCV(1) = OCVfromSOCtemp(z0,temp(1),model);
% 
%   % BEGIN MODIFYING CODE AFTER THIS
% 
%   T =temp(:); % sample code uses only single temperature -- you will need to change this!
%   
%   % The code reproduced below as a starting point for you is extracted from "simCell.m"
%   %
%   % Get model parameters from model structure -- notice that these retreive parameter values
%   % for only a single temperature. You will need to change this to load parameter values for
%   % all temperatures in temp.
%   
%   
%   RCfact = exp(-deltaT./abs(getParamESC('RCParam',T(1),model)))';
%   
%   if length(RCfact) ~= Nr,
%     error('iR0 does not have the correct number of entries');
%   end
%   G = zeros(size(T));Q = zeros(size(T));M = zeros(size(T));M = zeros(size(T)); M0 = zeros(size(T));
%   RParam = zeros(size(T));R0Param = zeros(size(T));etaParam = zeros(size(T));
%   for k=1:N,
%     G(k) = getParamESC('GParam',T(k),model); 
%     Q(k) = getParamESC('QParam',T(k),model);
%     etaParam(k) = getParamESC('etaParam',T(k),model);
%     M(k) = getParamESC('MParam',T(k),model);
%     M0(k) = getParamESC('M0Param',T(k),model);
%     RParam(k) = getParamESC('RParam',T(k),model);
%     R0Param(k) = getParamESC('R0Param',T(k),model);
%   end
%   etaik = ik; etaik(ik<0) = etaParam(ik<0).*ik(ik<0);
%   size(etaik)
%   size(G)
%   % Simulate the dynamic states of the model
%   for k = 2:length(ik),
%     rck(k,:) = rck(k-1,:)*diag(RCfact) + (1-RCfact')*etaik(k-1);
%      %rck(k,:) = rck(k-1,:)*RCfact(Nr,1)+ (1-RCfact(Nr,1))*etaik(k-1);
%   end
% 
%   for k =2:N,
%       zk(k) = zk(k-1)-etaik(k)*deltaT/(Q(k)*3600);
%   end
%   if any(zk>1.1),
%     warning('Current may have wrong sign as SOC > 110%');
%   end
%   
%   % Hysteresis stuff
%   for k=1:N
%     fac(k)=exp(-abs(G(k)*etaik(k)*deltaT/(3600*Q(k))));
%   end
%   
%   for k=2:length(ik),
%     hk(k)=fac(k-1)*hk(k-1)+(fac(k-1)-1)*sign(ik(k-1));
%     sik(k) = sign(ik(k));
%     if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
%   end
%     
%   % Compute output equation
%   for k= 1:N
%     OCV(k) = OCVfromSOCtemp(zk(k),T(k),model);
%   end
%   for k=1:N
%     vest(k,1) = OCV(k,1) - rck(k,1)*RParam(k,1) - R0Param(k,1)*ik(k,1) + M(k,1)*hk(k,1) + M0(k,1)*sik(k,1);
%   end
%   % FINISH MODIFYING CODE BEFORE THIS
% 
% end
% % END GRADED FUNCTION