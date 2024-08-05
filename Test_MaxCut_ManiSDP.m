function out = Test_MaxCut_ManiSDP(varargin)
    p = inputParser;
    addOptional(p, 'graph', 'G1', @ischar);
    addOptional(p, 'p0', 10, @isnumeric);
    addOptional(p, 'tol', 0.01, @isnumeric);
    addOptional(p, 'path2folder', '~/', @ischar);
    parse(p, varargin{:});

    %% Generate the max-cut problem
    maxcut_data = p.Results.graph;
    tol = p.Results.tol; % stopping tolerance
    path2folder = p.Results.path2folder; % path to the folder of this repo
    p0 = p.Results.p0; 
    data = load([path2folder, 'datasets/graphs/MaxCut/', maxcut_data, '.mat']);
    A = data.A;

    n = size(A,1);
    C = spdiags(A*ones(n,1),0,n,n) - A;
    C = 0.5*(C+C'); % symmetrize if not symmetric
    C = (-0.25).*C;

    %% Solve using ManiSDP
    rng(0);
    clear options;
    options.p0 = p0;
    options.tol = tol;
    options.gradtol = tol;
    options.delta = 4;
    tic
    [~, fval, data] = ManiSDP_onlyunitdiag(C, options);
    emani = data.dinf;
    tmani = toc;

    fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
end