function output = getInitParamsGalPathway( param )
% this function is used to get the initial value for ON/OFF state

load_global;
opt = odeset('NonNegative',1:12);
accurate_thresh = 10^-8;

% generate seed param for starting from GLUCOSE only state
param.exglu = 2*perc_to_nm;
param.exgal = 0;
odefunc = @(t,y)GALode5(t,y,param);

tmp = ones(1,12);
tmp(1) = 0;     % Gal1
tmp(6) = 0;     % Gal3*
tmp(8) = 0;     % C83
tmp(11) = 0;    % intracellular galactose
[~, y] = ode15s( odefunc, [0 10000], tmp, opt);
y(y<accurate_thresh) = 0;   % omit values that are too small
y0_Glu=y(end,:);

% generate seed param for starting from GALACTOSE only state
param.exglu = 0;
param.exgal = 2*perc_to_nm;
odefunc = @(t,y)GALode5(t,y,param);

tmp = ones(1,12);
tmp(1) = 0;     % Gal1
tmp(10) = 0;    % intracellular glucose
[~, y] = ode15s( odefunc, [0 10000], tmp, opt);
y(y<accurate_thresh) = 0;   
y0_Gal=y(end,:);
y0_Gal(1) = y0_Glu(1);

output.y0_Gal = max(y0_Gal,0);
output.y0_Glu = max(y0_Glu,0);

end

