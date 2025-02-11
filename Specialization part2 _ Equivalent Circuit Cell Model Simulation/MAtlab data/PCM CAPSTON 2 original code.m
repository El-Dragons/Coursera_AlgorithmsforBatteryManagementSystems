

end
% END GRADED FUNCTION
% GRADED FUNCTION (do not modify this line)

% function [vpack,vcell,icell,zcell,qcell,rcell] = simPCMTemp(Ns,Np,current,temp,deltaT,model)
%
% Simulate parallel-connected-module packs (cells are connected in parallel
% to make modules; these modules are connected in series to make packs).
% Note: This function must work for models having a single R-C pair. It is
% not required to work for models having multiple R-C pairs.
%
% Ns - number of modules connected in series to make a pack
% Np - number of cells connected in parallel in each module
% current - battery pack current, where (+) is discharge. Size is N x 1.
% temp  - temperature (degC). Size is N x 1 (all cells at same temperature).
% deltaT = sampling interval in data (s)
% model - standard model structure
% delta - variability vector: variable initial SOC if delta(1)==1; 
%                         variable total capacity if delta(2)==1;
%                         variable series resistance if delta(3)==1;
%
% vpack - battery pack voltage. Size is N x 1.
% vcell - individual cell voltages. Size is N x Ns x Np
% icell - individual cell currents. Size is N x Ns x Np
% zcell - individual cell states of charge. Size is N x Ns x Np
% qcell - individual cell capacities. Size is N x Ns x Np
% rcell - individual cell series resistances. Size is N x Ns x Np

function [vpack,vcell,icell,zcell,qcell,rcell] = simPCMTemp(Ns,Np,current,temp,deltaT,model,delta)

  % Force current to be column vector in case user entered data incorrectly
  current = current(:); N = length(current); temp = temp(:);
  % Initialize function outputs
  vpack = zeros(N,1); vcell = zeros(N,Ns,Np); icell = zeros(N,Ns,Np);
  zcell = zeros(N,Ns,Np); qcell = zeros(N,Ns,Np); rcell = zeros(N,Ns,Np);
  
  % Do some error checking on the function inputs
  Nr = length(getParamESC('RCParam',25,model)); % number of R-C pairs.
  if Nr ~= 1,
    error('This code does not work for models having multiple R-C pairs.');
  end
  if length(temp) ~= N,
    error('Input "temp" vector not the correct dimension.');
  end
  
  % Initialize states for ESC cell model
  if delta(1),
    z = reshape(linspace(0.3,0.7,Ns*Np),Ns,Np); % Different initial SOCs
  else
    z = 0.5*ones(Ns,Np);
  end
  irc = zeros(Ns,Np);
  h   = zeros(Ns,Np);

  % BEGIN MODIFYING CODE AFTER THIS
 
  T = temp(1); % sample code uses only single temperature -- you will need to change this!
  
  % The code reproduced below as a starting point for you is modified from
  % "simPCM.m".
  %
  % Get model parameters from model structure -- notice that these retreive
  % parameter values for only a single temperature. You will need to change
  % this to load parameter values for all temperatures in temp.
  
  % Default initialization for cells within the pack
  q  = getParamESC('QParam',T(1),model)*ones(Ns,Np); 
  rc = exp(-deltaT./abs(getParamESC('RCParam',T(1),model)))'*ones(Ns,Np);
  r  = (getParamESC('RParam',T(1),model))';
  m  = getParamESC('MParam',T(1),model)*ones(Ns,Np);
  g  = getParamESC('GParam',T(1),model)*ones(Ns,Np);
  r0 = getParamESC('R0Param',T(1),model)*ones(Ns,Np); 
  rt = 0.000125; % 125 microOhm resistance for each tab

 % for k=1:N,
%    q(k,:,:)  = getParamESC('QParam',T(k),model)*ones(Ns,Np); 
%    rc(k,:,:) = exp(-deltaT./abs(getParamESC('RCParam',T(k),model)))'*ones(Ns,Np);
%    r(k,:,:)  = (getParamESC('RParam',T(k),model))';
%    m(k,:,:)  = getParamESC('MParam',T(k),model)*ones(Ns,Np);
%    g(k,:,:)  = getParamESC('GParam',T(k),model)*ones(Ns,Np);
 %   r0(k,:,:) = getParamESC('R0Param',T(k),model)*ones(Ns,Np); 

%  end

  rt = 0.000125; % 125 microOhm resistance for each tab
  % How to modify capacity at any temperature if that input option is set... 
  % Don't change this functionality since the grader assumes it.

  if delta(2), 
    q = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*q; 
  end

  % How to modify resistances at any temperature if that input option is set...
  % Don't change this functionality since the grader assumes it.
  if delta(3),
    r0 = reshape(linspace(0.95,1.05,Ns*Np),Ns,Np).*r0; 
  end
  
  r0 = r0 + 2*rt; % add tab resistance to cell resistance
  % Okay... now to simulate pack performance using ESC cell model.
  
  for k = 1:N,
    v = OCVfromSOCtemp(z,T,model); % get OCV for each cell: Ns * Np matrix
    v = v + m.*h - r.*irc; % add in capacitor voltages and hysteresis
    V = (sum(v./r0,2) - current(k))./sum(1./r0,2);
    ik= (v-repmat(V,1,Np))./r0;                              
    z = z - (1/3600)*ik./q;  % Update each cell SOC
    irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
    fac = exp(-abs(g.*ik)./(3600*q));                        
    h = fac.*h + (fac-1).*sign(ik); % Update hysteresis voltages
    vpack(k)     = sum(V); % Store pack voltage                         
    vcell(k,:,:) = v - ik.*r0; % Store cell voltages
    zcell(k,:,:) = z; % Store cell SOCs
    icell(k,:,:) = ik; % Store cell currents
    qcell(k,:,:) = q; % Store cell capacities
    rcell(k,:,:) = r0; % Store cell resistances
  end % for k

  % FINISH MODIFYING CODE BEFORE THIS

end
% END GRADED FUNCTION