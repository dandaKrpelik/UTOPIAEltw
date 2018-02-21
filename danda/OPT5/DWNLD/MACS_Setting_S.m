%% MACS2 SETTINGS
np=6;
no=2;
name='DragCostUTOPIAE';
vlb0=zeros(1,np);
vub0=ones(1,np);

opt.popsize=500;         % popsize (forM4 each archive)
opt.maxnfeval=opt.popsize*ngmax;     % maximum number of f evals
opt.rhoini=1;           % initial span of each local hypercube (1=full domain)
opt.F=0.9;              % F, the parameter for Differential Evolution
opt.CR=0.9;             % CR, crossover probability
opt.p_social=0.2;         % popratio
opt.max_arch = 100;

if no==3
    
    opt.max_arch=150;          % archive size
    
end

opt.coord_ratio=1;
opt.contr_ratio=0.5;         % contraction ratio
opt.draw_flag=0;          % draw flag
opt.cp=0;          % constraints yes/no
opt.MBHflag=0;          % number of MBH steps
opt.explore_DE_strategy = 'rand';
opt.social_DE_strategy = 'DE/current-to-rand/1';
opt.v = 0;
opt.dyn_pat_search = 1;
opt.upd_subproblems = 0;
opt.max_rho_contr = 5;
opt.pat_search_strategy = 'standard';
opt.optimal_control = 0;
opt.vars_to_opt = ones(length(vlb0),1);
opt.func=[];
opt.cfunc=[];
