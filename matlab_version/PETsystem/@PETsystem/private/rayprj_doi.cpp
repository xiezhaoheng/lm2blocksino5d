#include "mex.h"
#include "rtracer.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs == 0)
        return;

#if USE_TOF
    mexErrMsgTxt("Not support TOF!");
#endif

#if USE_FPROJ
#if USE_BPROJ
    mexErrMsgTxt("USE_FPROJ and USE_BPROJ can't be together!");
    return;
#endif
#endif
    
#if USE_BPROJ
#if USE_FPROJ
    mexErrMsgTxt("USE_FPROJ and USE_BPROJ can't be together!");
    return;
#endif
#endif
    
#if USE_FPROJ

    // forward projector:
    //
    //  proj = rayfp(image, imgsize, voxsize, xtal_xy_position, ring_z_location, maxstep, lmdata, steps)
    //
    // xtal_xy_position: 2 x nXtals x nDOI
    // lmdata: [xtal1,ring1,doi1,xtal2,ring2,doi2,step]
    //

    // read in image content
    IMAGE_DATA_TYPE* image = 0;
    if (nrhs < 1) {
        mexErrMsgTxt("need an image!");
        return;
    } else {
        image = (IMAGE_DATA_TYPE*)mxGetPr(prhs[0]);
        if (image == NULL) {
            mexErrMsgTxt("invalid image!");
            return;
        }        
    }
    
    // size based on image content
    INT n0 = mxGetN(prhs[0]);
    INT n1 = mxGetM(prhs[0]);
    INT nv = n0 * n1;

#endif
    
#if USE_BPROJ
    
    REAL* proj = 0;
    if (nrhs < 1)  {
        mexErrMsgTxt("need projections!");
        return;
    } else {
        proj = mxGetPr(prhs[0]);
        if (proj == NULL) {
            mexWarnMsgTxt("no proj found!");
        }
    }

    INT num_of_projs = mxGetN(prhs[0]) * mxGetM(prhs[0]);
#endif
    
    // read in image dimension
    REAL* imgsize = 0; 
    if (nrhs < 2) {
        mexErrMsgTxt("not specify image dimension!");
        return;
    } else {
        imgsize = mxGetPr(prhs[1]);
        if ((imgsize == NULL) ||
            (imgsize[0] <= 0) || 
            (imgsize[1] <= 0) || 
            (imgsize[2] <= 0)) {
            mexErrMsgTxt("invalid image size!");
            return;
        }
    }     

    // convert to integers
    INT ni = INT(imgsize[0]);
    INT nj = INT(imgsize[1]);
    INT nk = INT(imgsize[2]);    
    INT nsize = ni * nj * nk;

#if USE_FPROJ    
    // compare to image vector size
    if (nsize != nv) {
        mexErrMsgTxt("invalid image vector or image size!");
        return;
    }
#endif

    // read in voxel size
    REAL* voxsize = 0;
    if (nrhs < 3) {
        mexErrMsgTxt("not specify voxel size!");    
        return;
    } else {
        voxsize = mxGetPr(prhs[2]);
        if ((voxsize == NULL) || 
            (voxsize[0] <= 0) || 
            (voxsize[1] <= 0) || 
            (voxsize[2] <= 0)) {
            mexErrMsgTxt("invalid voxel size!");
            return;
        }
    }
             
    // xtal transaxial position (2d) 
    XTAL_XY_POSITION* xytemp = 0;
    std::vector<std::vector<XTAL_XY_POSITION> > xtal_xy_position;
    if (nrhs < 4) {
        mexErrMsgTxt("not specify crystal transaxial coordinates!");
        return;
    } else {

        xytemp = (XTAL_XY_POSITION*)mxGetPr(prhs[3]);
        if ((xytemp == NULL) || 
            (mxGetM(prhs[3]) != 2)) {
            mexErrMsgTxt("invalid crystal positions!");
            return;      
        }
        
    }
    
    const INT* sz = mxGetDimensions(prhs[3]);
    INT num_of_xtals = (mxGetNumberOfDimensions(prhs[3]) == 3) ? sz[2] : sz[1];
    INT num_of_dois = (mxGetNumberOfDimensions(prhs[3]) == 3) ? sz[1] : 1;
    xtal_xy_position.resize(num_of_xtals);
    for (INT d = 0; d < num_of_xtals; d ++) {            
        xtal_xy_position[d].resize(num_of_dois);
        for (INT c = 0; c < num_of_dois; c ++) {
            
            xtal_xy_position[d][c].x = xytemp[c + d*num_of_dois].x;
            xtal_xy_position[d][c].y = xytemp[c + d*num_of_dois].y;
            
        }            
    }
    
    // ring z locations
    REAL* ring_z_locations = 0;
    if (nrhs < 5) {
        mexErrMsgTxt("not specify ring axial coordinates!");
        return;
    } else {
        ring_z_locations = mxGetPr(prhs[4]);
        if ((ring_z_locations == NULL)) {
            mexErrMsgTxt("invalid ring axial coordinates!");
            return;
        }
    }
    INT nring = mxGetM(prhs[4]) * mxGetN(prhs[4]);
    if (nring == 0) {
        mexErrMsgTxt("invalid ring axial coordinates!");
        return;
    }

    // read listmode data
    if (nrhs < 6) {
        mexErrMsgTxt("no listmode data found!");
        return;
    } else {
        
        if ((mxGetClassID(prhs[5]) != mxUINT16_CLASS) ||
            (mxGetPr(prhs[5]) == NULL)) {
            mexErrMsgTxt("invalid listmode data, must be UINT8!");
            return;
        }
        
        
        if (mxGetM(prhs[5]) != 6) {
            mexErrMsgTxt("invalid listmode format, must be 6 x num_of_prompts!");
            return;
        }
    }
    DOILMEVENT* lmdata = (DOILMEVENT*)mxGetPr(prhs[5]);
    INT num_of_prompts = mxGetN(prhs[5]);

#if USE_BPROJ
    
    if ((num_of_prompts != num_of_projs) && (proj != 0)) {
        mexErrMsgTxt("wrong number of projections or prompts!");
        return;
    }

#endif
    
    // show something
#if USE_FPROJ
    mexPrintf("forward projection ... ");
#endif
#if USE_BPROJ
    mexPrintf("backprojection ... ");
#endif
#if 1
    mexPrintf("\tImage size: %d,%d,%d\n"
              "\tVoxel size: %g,%g,%g\n"
              "\tNumber of crystals on a ring: %d\n"
              "\tNumber of DOI bins per crystals: %d\n"
              "\tNumber of rings: %d\n"
              "\tNumber of prompts: %d\n",
              ni, nj, nk, voxsize[0], voxsize[1], voxsize[2],
              num_of_xtals, num_of_dois, nring, num_of_prompts);
#endif

    // create raytracer
    ImageRayTracer raytracer(ni, nj, nk, 
                             voxsize[0], voxsize[1], voxsize[2]);


#if USE_FPROJ    
    // create projection buffer
    plhs[0] = mxCreateDoubleMatrix(num_of_prompts, 1, mxREAL);
    REAL* proj = mxGetPr(plhs[0]);
    
    INT n;
    DOILMEVENT e;
#if USE_OMP
    #pragma omp parallel for private(n, e)
#endif
    for (n = 0; n < num_of_prompts; n ++) {
        
        e = lmdata[n];
        
        REAL x1 = xtal_xy_position[e.xtal_id_1][e.doi_id_1].x;
        REAL y1 = xtal_xy_position[e.xtal_id_1][e.doi_id_1].y;
        REAL x2 = xtal_xy_position[e.xtal_id_2][e.doi_id_2].x;
        REAL y2 = xtal_xy_position[e.xtal_id_2][e.doi_id_2].y;
                    
        proj[n] = raytracer.fproj(x1, y1, ring_z_locations[e.ring_id_1],
                                  x2, y2, ring_z_locations[e.ring_id_2], 
                                  image);
    }
#endif
    
#if USE_BPROJ
    plhs[0] = mxCreateDoubleMatrix(nsize, 1, mxREAL);
    REAL* bpimg = mxGetPr(plhs[0]);
    
    INT n;
#if USE_OMP    
    #pragma omp parallel for private(n)
#endif    
    for (n = 0; n < num_of_prompts; n ++) {
        
        DOILMEVENT e = lmdata[n];
        REAL weight = (proj == 0) ? 1.0 : proj[n]; 
                
        if (weight > 0) {
            
            REAL x1 = xtal_xy_position[e.xtal_id_1][e.doi_id_1].x;
            REAL y1 = xtal_xy_position[e.xtal_id_1][e.doi_id_1].y;
            REAL x2 = xtal_xy_position[e.xtal_id_2][e.doi_id_2].x;
            REAL y2 = xtal_xy_position[e.xtal_id_2][e.doi_id_2].y;
            
            raytracer.bproj(x1, y1, ring_z_locations[e.ring_id_1],
                            x2, y2, ring_z_locations[e.ring_id_2], 
                            weight, bpimg);
        }
    }
    
#endif
    
    mexPrintf("OK, done!\n");

    return;    
}