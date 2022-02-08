#include "mex.h"
#include "rtracer.h"

/*
	Created by Jian Zhou
	Date: Sept. 2014
*/

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct  tagXTALPAIR {
    short id1;
    short id2;
}XTALPAIR;

typedef struct tagPOINT3D {
    REAL x;
    REAL y;
    REAL z;
    REAL w;
}POINT3D;

// KeV
// 100, 150, 200, 300, 400, 500, 600, 800,

static REAL ENERGY_MIN = 100.0;
static REAL ENERGY_MAX = 800.0;

// mm^-1
// 1.707E-02, 1.505E-02, 1.370E-02, 1.186E-02, 1.061E-02, 9.687E-03, 8.956E-03, 7.865E-03,
// water!
static REAL g_lut_energy_to_attn[] = {
0.0170700000000000, 0.0168132532528914, 0.0165720498118463, 0.0163448037931792, 0.0161301513673779, 
0.0159269119147616, 0.0157340572268371, 0.0155506868472584, 0.0153760081505274, 0.0152093201144653, 
0.0150500000000000, 0.0148896443837830, 0.0147360085288400, 0.0145886136836971, 0.0144470284188588, 
0.0143108627184954, 0.0141797629548658, 0.0140534075930964, 0.0139315035035928, 0.0138137827826245, 
0.0137000000000000, 0.0135801971550087, 0.0134642915572912, 0.0133520680675073, 0.0132433280772461, 
0.0131378878971601, 0.0130355773341223, 0.0129362384316644, 0.0128397243519247, 0.0127458983806274, 
0.0126546330393490, 0.0125658092916143, 0.0124793158312822, 0.0123950484432879, 0.0123129094281769, 
0.0122328070830114, 0.0121546552322152, 0.0120783728027554, 0.0120038834387755, 0.0119311151514049, 
0.0118600000000000, 0.0117843475306916, 0.0117103961183849, 0.0116380809900828, 0.0115673408175308, 
0.0114981174835976, 0.0114303558678168, 0.0113640036492513, 0.0112990111250480, 0.0112353310432223, 
0.0111729184483696, 0.0111117305391346, 0.0110517265363910, 0.0109928675611915, 0.0109351165216394, 
0.0108784380079208, 0.0108227981948075, 0.0107681647510084, 0.0107145067548075, 0.0106617946154763, 
0.0106100000000000, 0.0105563783130117, 0.0105036806132120, 0.0104518799987263, 0.0104009506629791, 
0.0103508678376940, 0.0103016077394905, 0.0102531475198122, 0.0102054652179421, 0.0101585397168828, 
0.0101123507018964, 0.0100668786215165, 0.0100221046508593, 0.0099780106570755, 0.0099345791667952, 
0.0098917933354329, 0.0098496369182260, 0.0098080942428943, 0.0097671501838110, 0.0097267901375906, 
0.0096870000000000, 0.0096456082694702, 0.0096047986064754, 0.0095645572384959, 0.0095248708485425, 
0.0094857265558832, 0.0094471118977600, 0.0094090148120360, 0.0093714236207177, 0.0093343270142996, 
0.0092977140368847, 0.0092615740720337, 0.0092258968293033, 0.0091906723314321, 0.0091558909021388, 
0.0091215431544982, 0.0090876199798620, 0.0090541125372963, 0.0090210122435052, 0.0089883107632159, 
0.0089560000000000, 0.0089225021979421, 0.0088894038503984, 0.0088566969778651, 0.0088243738234437, 
0.0087924268448955, 0.0087608487070415, 0.0087296322744883, 0.0086987706046654, 0.0086682569411569, 
0.0086380847073132, 0.0086082475001300, 0.0085787390843795, 0.0085495533869833, 0.0085206844916149, 
0.0084921266335199, 0.0084638741945450, 0.0084359216983646, 0.0084082638058968, 0.0083808953108993, 
0.0083538111357371, 0.0083270063273144, 0.0083004760531627, 0.0082742155976785, 0.0082482203585037, 
0.0082224858430422, 0.0081970076651070, 0.0081717815416917, 0.0081468032898615, 0.0081220688237575, 
0.0080975741517108, 0.0080733153734603, 0.0080492886774715, 0.0080254903383494, 0.0080019167143448, 
0.0079785642449476, 0.0079554294485641, 0.0079325089202760, 0.0079097993296760, 0.0078872974187783,
0.0078650000000000,
};

////////////////////////////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs == 0)
        return; 
    
    // single scatter simulator:
    //
    //  S = sss(emis_image, attn_map, imgsize, voxsize, xtal_xy_position, ring_z_location)
    //
    //
        
    // read in emiss_image content
    IMAGE_DATA_TYPE* emis_image = 0;
    if (nrhs < 1) {
        mexErrMsgTxt("need an emission image!");
        return;
    } else {
        emis_image = (IMAGE_DATA_TYPE*)mxGetPr(prhs[0]);
        if (emis_image == NULL) {
            mexErrMsgTxt("invalid emission image!");
            return;
        }        
    }
    
    IMAGE_DATA_TYPE* attn_map = 0;
    if (nrhs < 2) {
        mexErrMsgTxt("need an attenuation map!");
        return;
    } else {
        attn_map = (IMAGE_DATA_TYPE*)mxGetPr(prhs[1]);
        if (attn_map == NULL) {
            mexErrMsgTxt("invalid attenuation map!");
            return;
        }        
    }

    if (mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1])) {
        mexErrMsgTxt("size of emis != size of attn");
        return;
    }
    
    // size based on image content
    INT n0 = mxGetN(prhs[0]);
    INT n1 = mxGetM(prhs[0]);
    INT nv = n0 * n1;
    
    // read in image dimension
    REAL* imgsize = 0; 
    if (nrhs < 3) {
        mexErrMsgTxt("not specify image dimension!");
        return;
    } else {
        imgsize = mxGetPr(prhs[2]);
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
    
    if (nsize != nv) {
        mexErrMsgTxt("invalid image vector or image size!");
        return;
    }
    
    // read in voxel size
    REAL* voxsize = 0;
    if (nrhs < 4) {
        mexErrMsgTxt("not specify voxel size!");    
        return;
    } else {
        voxsize = mxGetPr(prhs[3]);
        if ((voxsize == NULL) || 
            (voxsize[0] <= 0) || 
            (voxsize[1] <= 0) || 
            (voxsize[2] <= 0)) {
            mexErrMsgTxt("invalid voxel size!");
            return;
        }
    }
    
    // xtal transaxial position (2d) 
    XTAL_XY_POSITION* xtal_xy_positions = 0;
    if (nrhs < 5) {
        mexErrMsgTxt("not specify crystal transaxial coordinates!");
        return;
    } else {
        xtal_xy_positions = (XTAL_XY_POSITION*)mxGetPr(prhs[4]);
        if ((xtal_xy_positions == NULL) || 
            (mxGetM(prhs[4]) != 2)) {
            mexErrMsgTxt("invalid crystal positions!");
            return;
        }
    }
    INT num_of_xtals = mxGetN(prhs[4]);
        
    // ring z locations
    REAL* ring_z_locations = 0;
    if (nrhs < 6) {
        mexErrMsgTxt("not specify ring axial coordinates!");
        return;
    } else {
        ring_z_locations = mxGetPr(prhs[5]);
        if ((ring_z_locations == NULL)) {
            mexErrMsgTxt("invalid ring axial coordinates!");
            return;
        }
    }
    INT nring = mxGetM(prhs[5]) * mxGetN(prhs[5]);
    if (nring == 0) {
        mexErrMsgTxt("invalid ring axial coordinates!");
        return;
    }
    
    // read in crystal pairs
    if (nrhs < 7) {
        mexErrMsgTxt("not specify xtal pairs!");
        return;
    } else {
        
        if ((mxGetClassID(prhs[6]) != mxINT16_CLASS) ||
            (mxGetPr(prhs[6]) == NULL)) {
            mexErrMsgTxt("invalid xtal pairs, must be INT16!");
            return;
        }
        
        if ((mxGetM(prhs[6]) != 2)) {
            mexErrMsgTxt("invalid xtal pairs!");
            return;
        }
    }
    INT num_of_xps = mxGetN(prhs[6]);
    XTALPAIR* xtal_pair = (XTALPAIR*)mxGetPr(prhs[6]);
        
    // read in TOF info
    REAL* TOF_info = 0;
    REAL tw_resolution = 0.0;
    REAL tw_spacing = 0.0;
    if (nrhs < 8) {
        mexErrMsgTxt("not specify TOF info!");
        return;
    } else {
        INT nd = mxGetM(prhs[7]) * mxGetN(prhs[7]);
        if (nd < 2) {
            mexErrMsgTxt("invalid TOF info, must supply timing resolution and timing bin size, both in ps!");
            return;
        } else {
            TOF_info = mxGetPr(prhs[7]);
            if (TOF_info == NULL) {
                mexErrMsgTxt("invalid TOF info, must supply timing resolution (FWHM) and timing bin size, both in ps!");
                return;
            }
            tw_resolution = TOF_info[0];
            tw_spacing = TOF_info[1];

            if (tw_resolution <= 0.0 ||
                tw_spacing <= 0.0) {
                mexErrMsgTxt("invalid TOF info!!");
                return;
            }
        }
    }

    int tbin_id = 0;
    if (nrhs < 9) {
        mexErrMsgTxt("not supply tbin_id! use 0");
        return;
    } else {
        if (prhs[8] == NULL) {
            mexErrMsgTxt("not supply tbin_id! use 0");
        } else {
            tbin_id = INT(mxGetScalar(prhs[8]));
        }
    }

#ifndef USE_TOF
    mexErrMsgTxt("USE_TOF must be on!");
    return;
#endif

// read energy window and resolution
    REAL* energy_info = 0;
    if (nrhs < 10) {
        mexErrMsgTxt("not specify energy info");
        return;
    } else {
        energy_info = mxGetPr(prhs[9]);
        if (energy_info == NULL) {
            mexErrMsgTxt("invalid energy info!");
            return;
        }

        if (mxGetM(prhs[9]) * mxGetN(prhs[9]) < 3) {
            mexErrMsgTxt("invalid energy info, must be [min, max, resolution] (unit in keV)");
            return;
        }
    }
    REAL energy_window_min = energy_info[0];
    REAL energy_window_max = energy_info[1];
    REAL energy_resolution = energy_info[2];
    REAL energy_sigma = energy_resolution / 2.3548;


#if 1
    mexPrintf("\tImage size: %d,%d,%d\n"
              "\tVoxel size: %g,%g,%g\n"
              "\tNumber of crystals on a ring: %d\n"
              "\tNumber of rings: %d\n"
              "\tNumber of LORs (per plane): %d\n"
              "\tTOF info: %f, %f (ps)\n"
              "\ttbin id: %d\n",
              ni, nj, nk, voxsize[0], voxsize[1], voxsize[2],
              num_of_xtals, nring, num_of_xps, tw_resolution, tw_spacing, tbin_id);
#endif

    // calc energy step size
    REAL energy_stepsize = (ENERGY_MAX - ENERGY_MIN) / 
                           (sizeof(g_lut_energy_to_attn) / sizeof(REAL)- 1);
      
    // create a raytracer    
    ImageRayTracer::TRUNCATION = 6;
    ImageRayTracer raytracer(ni, nj, nk, 
                             voxsize[0], voxsize[1], voxsize[2],
                             tw_resolution * 0.15 / (2*sqrt(2*log(2))),
                             tw_spacing * 0.15);

    REAL vk = voxsize[2];
    REAL vj = voxsize[1];
    REAL vi = voxsize[0];

    //
    // pick scatter points    
    // need to design
    //
    mexPrintf("Select scatter points ...\n");    
    REAL threshold = 1.0;
    std::vector<POINT3D> scatter_point_buffer;
    for (INT k = 0; k < nk; k ++) {
        REAL sz = (-nk*0.5 + 0.5 + k) * vk;
        for (INT j = 0; j < nj; j ++) {
            REAL sx = (-nj*0.5 + 0.5 + j) * vj;
            for (INT i = 0; i < ni; i ++) {
                REAL sy = (ni*0.5 - 0.5 - i) * vi;
                
                REAL prob = (rand() % RAND_MAX) / (REAL)RAND_MAX;
                
                if ((attn_map[i + j * ni + k * ni * nj] > 0.00001) &&
                    (prob < threshold)) {
                    
                    POINT3D pt;
                    pt.x = sx;
                    pt.y = sy;
                    pt.z = sz;
                    pt.w = attn_map[i + j * ni + k * ni * nj];
                    scatter_point_buffer.push_back(pt);
                                        
                }                
            }
        }
    }
    
    mexPrintf("# of selected scatter points: %d\n", scatter_point_buffer.size()); 
    
    //
    // precompute ray sum
    //
    mexPrintf("Precomputing ray-sum lookup table ...\n");
    std::vector<std::vector<std::vector<REAL> > > attn_ray_sum(nring);
    for (INT rid = 0; rid < nring; rid ++) {
        
        attn_ray_sum[rid].resize(num_of_xps);

        INT n = 0;
#if USE_OMP
        #pragma omp parallel for private(n)
#endif        
        for (n = 0; n < num_of_xtals; n ++) {
            
            REAL xtal_x = xtal_xy_positions[n].x;
            REAL xtal_y = xtal_xy_positions[n].y;
            REAL xtal_z = ring_z_locations[rid];
            
            attn_ray_sum[rid][n].resize(scatter_point_buffer.size());
            
            for (INT m = 0; m < scatter_point_buffer.size(); m ++) {
                
                REAL sx = scatter_point_buffer[m].x;                    
                REAL sy = scatter_point_buffer[m].y;
                REAL sz = scatter_point_buffer[m].z;
                
                attn_ray_sum[rid][n][m] = 
                    raytracer.fproj_nonTOF(sx, sy, sz, xtal_x, xtal_y, xtal_z, attn_map);

            }
            
        }
        
    }

    REAL ENERGY_CUT_MIN = 440;
    REAL ENERGY_CUT_MAX = 665;
         
    plhs[0] = mxCreateDoubleMatrix(num_of_xps, 1, mxREAL);
    REAL* S = mxGetPr(plhs[0]);
    
    for (INT r1 = 0; r1 < nring; r1 ++) {
    
        for (INT r2 = 0; r2 < nring; r2 ++) {

            INT n = 0;
#if USE_OMP
            #pragma omp parallel for private(n)
#endif            
            for (n = 0; n < num_of_xps; n ++) { 
                
//                if ((n+1) % 1000 == 0) {
                    //mexPrintf("processing #%d ...\n", n + 1);
                    //mexEvalString("drawnow;");
//                    printf("processing #%d ...\n", n + 1);
//                }
    
                XTALPAIR xp = xtal_pair[n];
                
                REAL xtal1_x = xtal_xy_positions[xp.id1-1].x;
                REAL xtal1_y = xtal_xy_positions[xp.id1-1].y;
                REAL xtal1_z = ring_z_locations[r1];

                REAL xtal2_x = xtal_xy_positions[xp.id2-1].x;
                REAL xtal2_y = xtal_xy_positions[xp.id2-1].y;
                REAL xtal2_z = ring_z_locations[r2];
                
                // distance between the two crystals
                REAL d12 = (xtal1_x - xtal2_x) * (xtal1_x - xtal2_x) + 
                           (xtal1_y - xtal2_y) * (xtal1_y - xtal2_y) + 
                           (xtal1_z - xtal2_z) * (xtal1_z - xtal2_z);                            
                
                S[n] = 0.0;
                                
                // loop over all scatter points   
                for (INT m = 0; m < scatter_point_buffer.size(); m ++) { //scatter_point_buffer.size()
                    
                    // scatter point (a.k. cone vertex)
                    REAL sx = scatter_point_buffer[m].x;                    
                    REAL sy = scatter_point_buffer[m].y;
                    REAL sz = scatter_point_buffer[m].z;                    
                    // attn at that point
                    REAL mu = scatter_point_buffer[m].w;
#if 0
mexPrintf("\nm=%d, A:[%f,%f,%f], S:[%f,%f,%f], B:[%f,%f,%f]\n", m, 
          xtal1_x, xtal1_y, xtal1_z, sx, sy, sz,
          xtal2_x, xtal2_y, xtal2_z);
#endif

                    // vector from xtal#1 to scatter point
                    REAL v1x = sx - xtal1_x;
                    REAL v1y = sy - xtal1_y;
                    REAL v1z = sz - xtal1_z;

                    // vector from scatter point to xtal#2
                    REAL v2x = xtal2_x - sx;
                    REAL v2y = xtal2_y - sy;
                    REAL v2z = xtal2_z - sz;
                    
                    // length of the two vectors
                    REAL len1_sq = v1x*v1x + v1y*v1y + v1z*v1z;
                    REAL len2_sq = v2x*v2x + v2y*v2y + v2z*v2z;
                            
                    // cosine of scattering angle
                    REAL cos_theta = (v1x * v2x + v1y* v2y + v1z * v2z) /
                                     (sqrt(len1_sq) * sqrt(len2_sq));
                       
                    // calc energy of the scattered photon
                    REAL E2 = 511.0 / (2 - cos_theta);

#if 0 // hard cut, no blurring                   
                    // apply energy cut (or energy resolution goes here)
                    if ((E2 > energy_window_max) || (E2 < energy_window_min)) {
                        continue;
                    }
                    REAL pb = 1.0;

#else // soft cut, compute probabilty within the energy window
                    REAL pb = 0.5 * (1.0 - erf((energy_window_min - E2)/(511.0 * energy_sigma))) - 
                              0.5 * (1.0 - erf((energy_window_max - E2)/(511.0 * energy_sigma)));
#endif

                    // ratio to determine attn coeff for scattered photon
                    // using water attn at various energies against the value at 511KeV
                    // linear interpolation based on lookup table
                    INT index = INT((E2 - ENERGY_MIN)/energy_stepsize);
                    REAL ratio = g_lut_energy_to_attn[index] / 0.0096; //0.009687;
                                    
                    // calc differential cross section, some constant coeffs are ignored
                    // using KN formula
                    REAL pn = 1.0 / (2 - cos_theta);
                    REAL KN = (pn*pn) * (pn + 1.0/pn - 1 + cos_theta*cos_theta) * 0.5;
                            
                    // ray intergral
                    REAL a1 = attn_ray_sum[r1][xp.id1-1][m];
                    REAL a2 = attn_ray_sum[r2][xp.id2-1][m];

//                    REAL e1 = emis_ray_sum[r1][xp.id1-1][m];
//                    REAL e2 = emis_ray_sum[r2][xp.id2-1][m];

                    REAL t = sqrt(len2_sq / len1_sq);
#if 0                    
mexPrintf("len:%f,%f, t=%f\n", sqrt(len1_sq), sqrt(len2_sq), t);
#endif
                    REAL vBx = sx + t * v1x;
                    REAL vBy = sy + t * v1y;
                    REAL vBz = sz + t * v1z;

#if 0
mexPrintf("[%f,%f,%f], virtual xtal2=[%f,%f,%f]\n", xtal1_x, xtal1_y, xtal1_z, vBx, vBy, vBz);
#endif

                    REAL e1 = raytracer.fproj_SSS_AS(xtal1_x, xtal1_y, xtal1_z,
                                                     sx, sy, sz, 
                                                     vBx, vBy, vBz, 
                                                     emis_image, tbin_id);

                    REAL inv_t = 1.0 / t;
                    REAL vAx = sx - inv_t * v2x;
                    REAL vAy = sy - inv_t * v2y;
                    REAL vAz = sz - inv_t * v2z;
#if 0
mexPrintf("virtual xtal1=[%f,%f,%f]\n", vAx, vAy, vAz);
#endif

                    REAL e2 = raytracer.fproj_SSS_SB(vAx, vAy, vAz,
                                                     sx, sy, sz,
                                                     xtal2_x, xtal2_y, xtal2_z, 
                                                     emis_image, tbin_id);

#if 0
mexPrintf("e1=%f,e2=%f\n\n", e1, e2);
#endif

                    REAL I = exp(-(a1 + ratio*a2)) * e1 + exp(-(a2 + ratio*a1)) * e2;                            
                            
                    // put together, ignore xtal efficiency (or assume identical to all)
                    // note: 1/(len1_sq*len2_sq) accounts for the solid angle effect
                    S[n] += 1.0 / (len1_sq*len2_sq) * mu * KN * I * pb;                         
                }
            }
        
        }
    }
    
    
    
}
