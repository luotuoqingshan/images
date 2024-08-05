function [outPGD,Xopt,yopt,Sopt] = Test_MaxCut_STRIDE(varargin)
    clc; %clear; close all; %restoredefaultpath;

    p = inputParser;
    addOptional(p, 'graph', 'G1', @ischar);
    addOptional(p, 'admm_tol', 0.1, @isnumeric);
    addOptional(p, 'tol', 0.01, @isnumeric);
    addOptional(p, 'path2folder', '~/', @ischar);
    parse(p, varargin{:});

    %% Generate the max-cut problem
    maxcut_data = p.Results.graph;
    admm_tol = p.Results.admm_tol; % tolerance for warmstart
    tol = p.Results.tol; % stopping tolerance
    path2folder = p.Results.path2folder; % path to the folder of this repo
    data = load([path2folder, 'datasets/graphs/MaxCut/', maxcut_data, '.mat']);
    
    %% Load problem data
    function SDP = maxcut_sdpt3(A)
        n = size(A, 1);
        SDP.blk = {'s', n};
        C = spdiags(A*ones(n,1),0,n,n) - A;
        C = 0.5*(C+C'); % symmetrize if not symmetric
        C = (-0.25).*C;
        SDP.C = {C};
        SDP.b = ones(n,1);
    
        rows = zeros(1, n);
        cols = zeros(1, n);
        vals = ones(1, n);
        for i = 1:n
            rows(1, i) = i * (i + 1)/2;
            cols(1, i) = i;
        end
        SDP.At = {sparse(rows, cols, vals, n*(n+1)/2, n)};
    end

    SDP = maxcut_sdpt3(data.A);
    
    

    %% Add STRIDE to matlab path, and provide path to dependencies
    addpath(genpath(pwd));
    manoptpath      = './manopt'; % required for local search
    sdpnalpath      = './SDPNAL+v1.0'; % required for ADMM+

    %% Solve using STRIDE
    addpath(genpath(manoptpath)); % add manopt to path

    % set parameters for STRIDE
    options.pgdStepSize     = 10; % step size, default 10
    options.maxiterPGD      = 5; % maximum outer iterations for STRIDE, default 5-10
    options.SDPNALpath      = sdpnalpath; % provide path to SDPNAL
    options.tolADMM         = admm_tol; % tolerance for warmstart, decrease this parameter for a better warmstart (but takes more time)
    options.tolPGD          = tol; % tolerance on KKT residual of the SDP
    options.lbfgseps        = false;

    % provide implementation to the local search method
    options.rrOpt           = 1:3; % round the leading 3 eigenvectors to generate hypotheses
    options.rrFunName       = 'local_search_quasar'; % name of the .m file that implements the local search

    % Primal initialization
    X0                  = [];

    % call STRIDE
    [outPGD,Xopt,yopt,Sopt] = PGDSDP(SDP.blk, SDP.At, SDP.b, SDP.C, X0, options);
    %infostride              = get_performance_quasar(Xopt,yopt,Sopt,SDP,R_gt);
end
