% make all 

% put number of threads
setenv('OMP_NUM_THREADS', '32');

% compile all
if 1
mex -o fproj_mt -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o fproj_mt_bresenham -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP -DUSE_BRESENHAM CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o fproj_mt_linterp -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP -DUSE_LINEAR_INTERP CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o fproj_tof_mt_bresenham -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP -DUSE_BRESENHAM -DUSE_TOF CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o fproj_tof_mt_linterp -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_FPROJ -DUSE_OMP -DUSE_TOF -DUSE_LINEAR_INTERP CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o bproj_mt -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_BPROJ -DUSE_OMP CXXOPTIMFLAGS="-O3 -DNDEBUG"


mex -o bproj_sq_mt -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_BPROJ -DUSE_OMP -DUSE_SQ_LOR CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o bproj_sq_tof_mt -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_BPROJ -DUSE_OMP -DUSE_SQ_LOR -DUSE_TOF CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o fproj_tof_mt -cxx -v rayprj.cpp -DUSE_FPROJ -DUSE_TOF -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o bproj_tof_mt -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_TOF -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3 -DNDEBUG"

mex -o bproj_mt_bresenham -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_BPROJ -DUSE_OMP -DUSE_BRESENHAM CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o bproj_mt_linterp -cxx -v rayprj.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_BPROJ -DUSE_OMP -DUSE_LINEAR_INTERP CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o bproj_tof_mt_bresenham -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_TOF -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_BRESENHAM CXXOPTIMFLAGS="-O3 -DNDEBUG"
mex -o bproj_tof_mt_linterp -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_TOF -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DUSE_LINEAR_INTERP CXXOPTIMFLAGS="-O3 -DNDEBUG"

end

if 0
mex -o fproj_st -cxx -v rayprj.cpp -DUSE_FPROJ 
mex -o bproj_st -cxx -v rayprj.cpp -DUSE_BPROJ 
mex -o bproj_st_bresenham -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_OMP -DUSE_BRESENHAM 
mex -o bproj_st_linterp -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_OMP -DUSE_LINEAR_INTERP 
mex -o fproj_tof_st -cxx -v rayprj.cpp -DUSE_FPROJ -DUSE_TOF -DCFOV_ENABLED 
mex -o bproj_tof_st -cxx -v rayprj.cpp -DUSE_BPROJ -DUSE_TOF -DCFOV_ENABLED 
%mex -o fproj_pawpet -cxx -v rayprj_pawpet.cpp -DUSE_FPROJ -DCFOV_ENABLED -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" 
%mex -o bproj_pawpet -cxx -v rayprj_pawpet.cpp -DUSE_BPROJ -DCFOV_ENABLED -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" 
%mex -o bproj_pawpet_rect -cxx -v rayprj_pawpet.cpp -DUSE_BPROJ -DUSE_OMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" 
end

