<<<<<<< HEAD
%% run refinement
function testGENFire_template(varargin)
%NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd
%ref_angles_NiPt_Mo_FS_171118_t2
optional_arg = {'-pr','data\small_projections.mat' ,...
                '-an','data\small_angles.mat' ,...
       '-r','3', '-m','1', '-d','0.5', '-it','100', '-ds','0.5', '-sm','0' };
optional_arg1 = optional_arg(1:2:end);
optional_arg2 = optional_arg(2:2:end);

if mod(nargin,2)~=0
    fprintf('not enough argument\n');
else
    varargin1 = varargin(1:2:nargin);
    varargin2 = varargin(2:2:nargin);
    
    for i=1:nargin/2
        arg_i = varargin1(i);
        [arg_exist,arg_ind] = ismember(arg_i,optional_arg1);
        if arg_exist, optional_arg2(arg_ind) = varargin2(i); end        
    end
end

pj_filename     = optional_arg2{1};
angle_filename  = optional_arg2{2};
oversamplingRatio   = eval(optional_arg2{3});
griddingMethod      = eval(optional_arg2{4});
distance        = eval(optional_arg2{5});
numIterations   = eval(optional_arg2{6});
ds              = eval(optional_arg2{7});
smooth          = eval(optional_arg2{8});

prefixstr = 'NiPt_Mo_181118BA_t2';
addpath('functions')
%/u/project/miao/minhrose/tomography/GENFIRE_matlab/
% absolute path on hoffman2

%pj_filename              = 'data/small_projections.mat';
%angle_filename           = 'data/small_angles.mat';

%pj_filename              = load('data/NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd.mat');
%angle_filename           = 'data\ref_angles_NiPt_Mo_FS_171118_t2.mat';

results_filename         = 'data/result.mat';
doGPU = 0;

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                         %%
%%%                        Welcome to GENFIRE!                              %%
%%%           GENeralized Fourier Iterative REconstruction                  %%
%%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Author: Alan (AJ) Pryor, Jr.
%%% email:  apryor6@gmail.com
%%% Jianwei (John) Miao Coherent Imaging Group
%%% University of California, Los Angeles
%%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isdeployed
    %addpath ../src/
    % addpath ../data/
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PLOT_YN = 0;

GENFIRE = GENFIRE_Reconstructor();

%%% See the README for description of parameters

GENFIRE.filename_Projections = pj_filename;
GENFIRE.filename_Angles = angle_filename ;
GENFIRE.filename_Results = results_filename;

GENFIRE.numIterations               = numIterations;    % num_iterations
GENFIRE.pixelSize = .5;
GENFIRE.oversamplingRatio           = oversamplingRatio;% ratio = 3
GENFIRE.griddingMethod              = griddingMethod;   % griddingMethod=2 for DFT
GENFIRE.allowMultipleGridMatches    = 1;
GENFIRE.constraintEnforcementMode   = 3;                      % no suppression
GENFIRE.interpolationCutoffDistance = distance;         % Cut-off distance
GENFIRE.dt_type                     = 1;                      % dt
GENFIRE.ds                          = ds;               % ds
GENFIRE.smooth                      = smooth;           % smoothness

GENFIRE.constraintPositivity = 1;
GENFIRE.constraintSupport = 1;
GENFIRE.ComputeFourierShellCorrelation = 0;
GENFIRE.numBins = 50;
GENFIRE.percentValuesForRfree = 0.05;
GENFIRE.numBinsRfree = 35;
GENFIRE.doCTFcorrection = 0;
GENFIRE.CTFThrowOutThreshhold = 0;
GENFIRE.calculate_Rfree = 0;
GENFIRE.DFT_doGPU = doGPU;


%GENFIRE.particleWindowSize = [];
%GENFIRE.phaseErrorSigmaTolerance = [];
%GENFIRE.vector1 = [0 0 1];
%GENFIRE.vector2 = [0 1 0];
GENFIRE.vector3 = [1 0 0];
GENFIRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin GENFIRE

% Read data from files. The data can also be directly be fed into GENFIRE
% Class, not necessarily from file
GENFIRE = readFiles(GENFIRE);

% Check if the data is consistent, and prepare reconstruction variables
% based on given parameters
GENFIRE = CheckPrepareData(GENFIRE);

% Run GENFIRE gridding
GENFIRE = runGridding(GENFIRE);

% Run FSC if flag is on
if GENFIRE.ComputeFourierShellCorrelation
    
    % Do not calculate R_free for FSC calculation
    calculate_Rfree_ori = GENFIRE.calculate_Rfree;
    GENFIRE.calculate_Rfree= 0;
    
    GENFIRE = runFSC(GENFIRE);
    
    % set the R_free flag back go original
    GENFIRE.calculate_Rfree=calculate_Rfree_ori;
end

GENFIRE = ClearCalcVariables(GENFIRE);

%run reconstruction
%GENFIRE = reconstruct(GENFIRE);
GENFIRE = reconstruct_dr(GENFIRE);
reconstruction = GENFIRE.reconstruction;
%final_rec = GENFIRE.final_rec;

%%
figure;img(permute(reconstruction,[2,3,1]));
figure;img(permute(reconstruction,[1,3,2]));
figure;img(reconstruction);
%%
%{
final_rec = GENFIRE.final_rec;
figure;img(permute(final_rec,[2,3,1]));
figure;img(permute(final_rec,[1,3,2]));
figure;img(final_rec);
%%
%save('Ni_DR.mat','final_rec');
%%
load('data\Up_NiPt_1118t2_FS.mat')
%%
%my_rec = tom_mrcread('results_pt_DR_07.mrc');
[SumX_1118t2_dr_dt01,SumY_1118t2_dr_dt01,SumZ_1118t2_dr_dt01,ESTvol_dr_dt01,peakX_1118t2_dr_dt01,peakY_1118t2_dr_dt01,peakZ_1118t2_dr_dt01,~] = Rot_FinalVol(reconstruction,Up);
figure;img(ESTvol_dr_dt01)
%}
=======
%% run refinement
function testGENFire_template(varargin)
%NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd
%ref_angles_NiPt_Mo_FS_171118_t2
optional_arg = {'-pr','/u/project/miao/minhrose/tomography/GENFIRE_matlab/data/NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd.mat' ,...
                '-an','/u/project/miao/minhrose/tomography/GENFIRE_matlab/data/ref_angles_NiPt_Mo_FS_171118_t2.mat' ,...
       '-r','4', '-m','2', '-d','0.1', '-it','200', '-ds','0.5', '-dt','1', '-sm','0' };
optional_arg1 = optional_arg(1:2:end);
optional_arg2 = optional_arg(2:2:end);

if mod(nargin,2)~=0
    fprintf('not enough argument\n');
else
    varargin1 = varargin(1:2:nargin);
    varargin2 = varargin(2:2:nargin);
    
    for i=1:nargin/2
        arg_i = varargin1(i);
        [arg_exist,arg_ind] = ismember(arg_i,optional_arg1);
        if arg_exist, optional_arg2(arg_ind) = varargin2(i); end        
    end
end

pj_filename         = optional_arg2{1};
angle_filename      = optional_arg2{2};
oversamplingRatio   = eval(optional_arg2{3});
griddingMethod      = eval(optional_arg2{4});
distance            = eval(optional_arg2{5});
numIterations       = eval(optional_arg2{6});
ds                  = eval(optional_arg2{7});
dt_type             = eval(optional_arg2{8});
smooth              = eval(optional_arg2{9});

prefixstr = 'NiPt_Mo_181118BA_t2';
%addpath('functions')
%/u/project/miao/minhrose/tomography/GENFIRE_matlab/
% absolute path on hoffman2

%pj_filename              = 'data/small_projections.mat';
%angle_filename           = 'data/small_angles.mat';

%pj_filename              = load('data/NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd.mat');
%angle_filename           = 'data\ref_angles_NiPt_Mo_FS_171118_t2.mat';

results_filename         = 'data/result.mat';
doGPU = 1;

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                         %%
%%%                        Welcome to GENFIRE!                              %%
%%%           GENeralized Fourier Iterative REconstruction                  %%
%%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Author: Alan (AJ) Pryor, Jr.
%%% email:  apryor6@gmail.com
%%% Jianwei (John) Miao Coherent Imaging Group
%%% University of California, Los Angeles
%%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isdeployed
    %addpath ../src/
    % addpath ../data/
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PLOT_YN = 0;

GENFIRE = GENFIRE_Reconstructor();

%%% See the README for description of parameters

GENFIRE.filename_Projections = pj_filename;
GENFIRE.filename_Angles = angle_filename ;
GENFIRE.filename_Results = results_filename;

GENFIRE.numIterations               = (numIterations);    % num_iterations
GENFIRE.pixelSize = .5;
GENFIRE.oversamplingRatio           = (oversamplingRatio);% ratio = 3
GENFIRE.griddingMethod              = (griddingMethod);   % griddingMethod=2 for DFT
GENFIRE.allowMultipleGridMatches    = 1;
GENFIRE.constraintEnforcementMode   = 3;                      % no suppression
GENFIRE.interpolationCutoffDistance = (distance);         % Cut-off distance
GENFIRE.dt_type                     = dt_type;                      % dt
GENFIRE.ds                          = (ds);               % ds
GENFIRE.smooth                      = (smooth);           % smoothness

GENFIRE.constraintPositivity = 1;
GENFIRE.constraintSupport = 1;
GENFIRE.ComputeFourierShellCorrelation = 0;
GENFIRE.numBins = 50;
GENFIRE.percentValuesForRfree = 0.05;
GENFIRE.numBinsRfree = 35;
GENFIRE.doCTFcorrection = 0;
GENFIRE.CTFThrowOutThreshhold = 0;
GENFIRE.calculate_Rfree = 0;
GENFIRE.DFT_doGPU = doGPU;


%GENFIRE.particleWindowSize = [];
%GENFIRE.phaseErrorSigmaTolerance = [];
%GENFIRE.vector1 = [0 0 1];
%GENFIRE.vector2 = [0 1 0];
GENFIRE.vector3 = [1 0 0];
GENFIRE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin GENFIRE

% Read data from files. The data can also be directly be fed into GENFIRE
% Class, not necessarily from file
GENFIRE = readFiles(GENFIRE);

% Check if the data is consistent, and prepare reconstruction variables
% based on given parameters
GENFIRE = CheckPrepareData(GENFIRE);

% Run GENFIRE gridding
GENFIRE = runGridding(GENFIRE);

% Run FSC if flag is on
if GENFIRE.ComputeFourierShellCorrelation
    
    % Do not calculate R_free for FSC calculation
    calculate_Rfree_ori = GENFIRE.calculate_Rfree;
    GENFIRE.calculate_Rfree= 0;
    
    GENFIRE = runFSC(GENFIRE);
    
    % set the R_free flag back go original
    GENFIRE.calculate_Rfree=calculate_Rfree_ori;
end

GENFIRE = ClearCalcVariables(GENFIRE);

%run reconstruction
%GENFIRE = reconstruct(GENFIRE);
GENFIRE = reconstruct_dr(GENFIRE);
reconstruction = GENFIRE.reconstruction;
%final_rec = GENFIRE.final_rec;
dir = pwd();
[~,job_ID] = system('echo $JOB_ID');
job_ID = job_ID(~isspace(job_ID));
filename = strcat(dir, '/pt_result_DR_', num2str(1/smooth), '_', job_ID, '.mat');

save(filename,'reconstruction');
%%
%{
figure;img(permute(reconstruction,[2,3,1]));
figure;img(permute(reconstruction,[1,3,2]));
figure;img(reconstruction);
%%
final_rec = GENFIRE.final_rec;
figure;img(permute(final_rec,[2,3,1]));
figure;img(permute(final_rec,[1,3,2]));
figure;img(final_rec);
%%
%save('Ni_DR.mat','final_rec');
%%
load('data\Up_NiPt_1118t2_FS.mat')
%%
%my_rec = tom_mrcread('results_pt_DR_07.mrc');
[SumX_1118t2_dr_dt01,SumY_1118t2_dr_dt01,SumZ_1118t2_dr_dt01,ESTvol_dr_dt01,peakX_1118t2_dr_dt01,peakY_1118t2_dr_dt01,peakZ_1118t2_dr_dt01,~] = Rot_FinalVol(reconstruction,Up);
figure;img(ESTvol_dr_dt01)
%}
>>>>>>> 633f82bdb29f44ac494198ed3ac17fc7d1999ce9
