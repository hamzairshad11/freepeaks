function build_ofec_mpmmo_mex()
%BUILD_OFEC_MPMMO_MEX Build the MATLAB wrapper for the OFEC MPMMO DLL.

root = fileparts(fileparts(mfilename('fullpath')));
matlabDir = fullfile(root,'Matlab');
dllDir = fullfile(root,'FreePeaks_Benchmark','FreePeaks_Benchmark','x64','Debug');
cmakeDllPath = fullfile(root,'build','ofec_mpmmo','bin','ofec_mpmmo.dll');
legacyDllPath = fullfile(dllDir,'ofec_mpmmo.dll');
if isfile(cmakeDllPath)
    dllPath = cmakeDllPath;
else
    dllPath = legacyDllPath;
end

if ~isfile(dllPath)
    error('OFEC:MPMMO:MissingDLL','Build ofec_mpmmo.dll first: %s',dllPath);
end

copyfile(dllPath,fullfile(matlabDir,'ofec_mpmmo.dll'));
mex('-v','-outdir',matlabDir,fullfile(matlabDir,'ofec_mpmmo_mex.cpp'));
fprintf('Built %s\n',fullfile(matlabDir,['ofec_mpmmo_mex.',mexext]));
end
