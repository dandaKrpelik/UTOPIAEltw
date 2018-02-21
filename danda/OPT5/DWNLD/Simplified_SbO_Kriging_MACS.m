clear all
close all
clc

addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/smart-o2c/Optimisation');
addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/smart-o2c/Optimisation/MACS');
addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/DWNLD');

rng('default')
rng(220)

disp("DANDA WAS HERE");

ndim=6;
%% DOE
nn=100;

xKept=lhsdesign(nn,ndim,'criterion','maximin','iterations',30);
TRUE_F=@(x) mask_DragCostUTOPIAE(x);   				%######## RESCALE and GIVES BACK obj func?
%% Evaluation of the two objectives
yKept1=[];
yKept2=[];
for i=1:nn
    [f]=TRUE_F(xKept(i,:));
    yKept1=[yKept1; f(1)];
    yKept2=[yKept2; f(2)];
end


%% Creation of the first Kriging models (one for each objective function)
addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/DWNLD/dace')
theta = [10 10 10 10 10 10]; lob = [1e-1 1e-1 1e-1 1e-1 1e-1 1e-1]; upb = [25 25 25 25 25 25];
[dmodel1, perf1] = dacefit(xKept,yKept1, @regpoly0, @corrgauss, theta, lob, upb);
[dmodel2, perf2] = dacefit(xKept,yKept2, @regpoly0, @corrgauss, theta, lob, upb);
dmodel10=dmodel1;
dmodel20=dmodel2;

%% Definition of the surrogate model that should be passed to optimiser
addpath('/home/utopiae_esr/Documents/UTOPIAEltw/danda/OPT5/DWNLD');
FITNESSFCN=@mask_surrogateK;		
ng=0;
ngmax=0;
while size(yKept1,1)<600 && ng<200
    %% Setting of the GA options
    PopSize=5000;
    X0=lhsdesign(PopSize,ndim,'criterion','maximin','iterations',30);
    memoryI=[];
    naddMax=20;
    naddZ=0;
    ngmax=ngmax+4;
    %% Surrogate based Optimisation loop
    while naddZ<=2 && size(yKept1,1)<600
        ng=ng+1;
        MACS_Setting_S
        
        disp('START MACS')
        [X,FVAL,exitflag,output] = optimise_macs(FITNESSFCN,memoryI,vlb0,vub0,opt,dmodel1,dmodel2);
        disp('END MACS')
        X0=X;
        POPULATION=X;
        
        update_database
        
        if nadd>0
            [dmodel1, perf1] = dacefit(xKept,yKept1, @regpoly0, @corrgauss, theta, lob, upb);
            [dmodel2, perf2] = dacefit(xKept,yKept2, @regpoly0, @corrgauss, theta, lob, upb);
            
            
            FITNESSFCN=@mask_surrogateK;

            Nfopt=[]
            for i=1:size(X0,1)
                Nfopt=[Nfopt;FITNESSFCN(X0(i,:),dmodel1,dmodel2)];
            end
            memoryI=[X0 Nfopt zeros(size(X0,1),2)];
        end
        if nadd==0
            naddZ=naddZ+1;
            memoryI=[X0 FVAL(:,2:3) zeros(size(X0,1),2)];
        else
            naddZ=0;
        end
        disp([size(yKept1) nadd naddZ])
        pause(2)
    end
end
