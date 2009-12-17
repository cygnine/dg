function[dg] = init__()
% init__ -- Initialization file for dg package
%
% [dg] = init__()

dg = recurse_files(pwd);
dg.fluxes = matlab_import('fluxes');
dg.meshes = matlab_import('meshes');
