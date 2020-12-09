function G = getG(param,method)
%%GETG
% w.r.t. the choice of method constructs dynamic system matrices ABCD 
% possible method choices: 
%   method =
%       'gd_pop'    : GD popular
%       'gd_opt'    : GD optimal
%       'nag_std'   : NAG standard NAG  : Nesterov Accelerated GD
%       'nag_opt'   : NAG optimal
%       'hbm_opt'   : HBM optimal       : Heavy Ball Method
%
% param structure should include the following fields
%   param.alpha
%   param.beta
%
%   G = getG(param,method)
%       where G.A, G.B, G.C. G.D
%
switch method
    %GD popular choice
    case 'gd_pop'
        G = G_gd(param.alpha);
        
    %NAG standard choice
    case 'nag_std' 
        G = G_nag(param.alpha,param.beta);
        
    %GD optimal tuning
    case 'gd_opt' 
        G = G_gd(param.alpha);
        
    %NAG optimal tuning
    case 'nag_opt'
        G = G_nag(param.alpha,param.beta);
        
    %HBM optimal tuning
    case 'hbm_opt' 
        G = G_hbm(param.alpha,param.beta);
end
end