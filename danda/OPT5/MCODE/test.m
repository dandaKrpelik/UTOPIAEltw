addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/smart-o2c/Optimisation/MACS');
addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/DWNLD');


no_trials = 0;
no_par = 6;
no_obj = 2;
% 
% SHEET = zeros(no_trials, no_par + no_obj)
% for i=1:no_trials
%    x = [ rand()*10
%           rand()*10
%            rand()*100
%            rand()*1000
%            rand()*2000
%            rand()*10000 ];
%    y = DragCostUTOPIAE(x);
%     
%    SHEET[i,1:no_par] = x[:];
%    SHEET[i,no_par+1:end] = y[:];
% end

Simplified_SbO_Kriging_MACS
