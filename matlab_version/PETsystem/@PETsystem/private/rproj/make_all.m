% make all 

% put number of threads
setenv('OMP_NUM_THREADS', '32');

% compile all
if 1
mex -o sss -cxx -v sss.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o sss_tof -cxx -v sss_tof.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP -DUSE_TOF CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o sss_tof_fast -cxx -v sss_tof_fast.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP -DUSE_TOF CXXOPTIMFLAGS="-O3 -DNDEBUG"
end

if 0
mex -o fproj_st -cxx -v rayprj.cpp -DUSE_FPROJ
mex -o bproj_st -cxx -v rayprj.cpp -DUSE_BPROJ
mex -o fproj_tof_st -cxx -v rayprj.cpp -DUSE_FPROJ -DUSE_TOF -DCFOV_ENABLED
mex -o bproj_tof_st -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_TOF -DCFOV_ENABLED
end

