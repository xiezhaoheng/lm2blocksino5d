#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

#ifndef RTRACER_H
#define RTRACER_H

/*
	Created by Jian Zhou
	Date: Jan 2012
*/

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
typedef int     INT;
typedef double  REAL;
typedef double  IMAGE_DATA_TYPE;

typedef struct  tagXTAL_XY_POSITION {
    double x;
    double y;
}XTAL_XY_POSITION;

typedef struct  tagLMEVENT {
    short xtal_id_1;
    short ring_id_1;
    short xtal_id_2;
    short ring_id_2;
    short tbin_id;
}LMEVENT;

typedef struct tagDOILMEVENT {
    unsigned char xtal_id_1;
    unsigned char ring_id_1;
    unsigned char doi_id_1;
    unsigned char xtal_id_2;
    unsigned char ring_id_2;
    unsigned char doi_id_2;
    unsigned char step_id;
}DOILMEVENT;

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ImageRayTracer
//
////////////////////////////////////////////////////////////////////////////////////////////////////
#define DIST_NO_SQRT
//#define BP_ATOMIC 1

class ImageRayTracer
{
public:
    ImageRayTracer(const INT vdim_i = 0,
                   const INT vdim_j = 0,
                   const INT vdim_k = 0,
                   const REAL vox_size_i = 0,
                   const REAL vox_size_j = 0,
                   const REAL vox_size_k = 0,
                   const REAL tw_sigma = -1,
                   const REAL tw_spacing = -1);
    ~ImageRayTracer();

public:
    REAL fproj(REAL p0x, REAL p0y, REAL p0z,
               REAL p1x, REAL p1y, REAL p1z,
               const IMAGE_DATA_TYPE* img, const INT tbin_id = 0);

    REAL fproj_SSS_AS(REAL, REAL, REAL,
                   REAL, REAL, REAL,
                   REAL, REAL, REAL,
                   const IMAGE_DATA_TYPE* img, const INT tbin_id);

    void fproj_TOF_SSS_AS(REAL p0x, REAL p0y, REAL p0z, // A
                          REAL sx, REAL sy, REAL sz, // S
                          REAL virtual_x, REAL virtual_y, REAL virtual_z, // B'
                          const IMAGE_DATA_TYPE* img, 
                          const std::vector<INT>& tbin, REAL* out);


    REAL fproj_SSS_SB(REAL, REAL, REAL,
                   REAL, REAL, REAL,
                   REAL, REAL, REAL,
                   const IMAGE_DATA_TYPE* img, const INT tbin_id);

    void fproj_TOF_SSS_SB(REAL virtual_x, REAL virtual_y, REAL virtual_z, // A
                          REAL sx, REAL sy, REAL sz, // S
                          REAL p1x, REAL p1y, REAL p1z, // B'
                          const IMAGE_DATA_TYPE* img, 
                          const std::vector<INT>& tbin, REAL* out);


    REAL fproj_nonTOF(REAL p0x, REAL p0y, REAL p0z,
                      REAL p1x, REAL p1y, REAL p1z,
                      const IMAGE_DATA_TYPE* img);
    
    void bproj(REAL p0x, REAL p0y, REAL p0z,
               REAL p1x, REAL p1y, REAL p1z,
               const REAL weight, IMAGE_DATA_TYPE* img, const INT tbin_id = 0);
        
private:
    bool hitCheck(const REAL p0x, const REAL p0y, const REAL p0z,
                  const REAL p1x, const REAL p1y, const REAL p1z,
                  REAL& t_min, REAL& t_max);
    void calcTOFLOREndPoints(const INT tbin_id,
                             REAL& p0x, REAL& p0y, REAL& p0z,
                             REAL& p1x, REAL& p1y, REAL& p1z,
                             REAL& nrm_dx, REAL& nrm_dy, REAL& nrm_dz);

    void calcTOFBinCenter(REAL p0x, REAL p0y, REAL p0z,
                          REAL p1x, REAL p1y, REAL p1z,
                          REAL& nrm_dx, REAL& nrm_dy, REAL& nrm_dz,
                          const std::vector<INT>& tbin,
                          std::vector<REAL>& tbin_center_list_x,
                          std::vector<REAL>& tbin_center_list_y,
                          std::vector<REAL>& tbin_center_list_z);

    
private:
    INT m_vdim_i;
    INT m_vdim_j;
    INT m_vdim_k;
    INT m_vdim_ixj;
    REAL m_vox_size_i;
    REAL m_vox_size_j;
    REAL m_vox_size_k;
    REAL m_vbd_x0;
    REAL m_vbd_x1;
    REAL m_vbd_y0;
    REAL m_vbd_y1;
    REAL m_vbd_z0;
    REAL m_vbd_z1;
    REAL m_tw_sigma;
    REAL m_tw_spacing;
    std::vector<REAL> m_gw_lut;
    REAL m_gs_inv;
    
public:
    static const REAL PI;
    static const REAL FLTMIN;
    static INT TRUNCATION;
    
private:
    static const INT GWSAMPLESIZE;
};

const REAL ImageRayTracer::FLTMIN = 1e-12;
const REAL ImageRayTracer::PI = 3.1415926535897932384626;
const INT ImageRayTracer::GWSAMPLESIZE = 5120;
INT ImageRayTracer::TRUNCATION = 3;

ImageRayTracer::ImageRayTracer(const INT vdim_i,
                               const INT vdim_j,
                               const INT vdim_k,
                               const REAL vox_size_i,
                               const REAL vox_size_j,
                               const REAL vox_size_k,
                               const REAL tw_sigma,
                               const REAL tw_spacing) :
    m_vdim_i(vdim_i),
    m_vdim_j(vdim_j),
    m_vdim_k(vdim_k),
    m_vox_size_i(vox_size_i),
    m_vox_size_j(vox_size_j),
    m_vox_size_k(vox_size_k),
    m_vdim_ixj(vdim_i * vdim_j),
    m_tw_sigma(tw_sigma),
    m_tw_spacing(tw_spacing)
{
    m_vbd_y0 = -(vdim_i * vox_size_i) * 0.5;
    m_vbd_y1 = +(vdim_i * vox_size_i) * 0.5;
    m_vbd_x0 = -(vdim_j * vox_size_j) * 0.5;
    m_vbd_x1 = +(vdim_j * vox_size_j) * 0.5;
    m_vbd_z0 = -(vdim_k * vox_size_k) * 0.5;
    m_vbd_z1 = +(vdim_k * vox_size_k) * 0.5;

#if USE_TOF    
    if (m_tw_sigma < 0 || m_tw_spacing < 0) {
        mexErrMsgTxt("invalid timing info!");
        abort();
    }

    // precompute gaussian window lookuptable
    // +/- 3-sigma truncation
#if 0  
    REAL gw_c0 = 1.0 / (2.0 * m_tw_sigma * m_tw_sigma);
    REAL gw_c1 = 1.0 / (sqrt(2.0 * PI) * m_tw_sigma);
    REAL gw_stepsize = 3.0 * m_tw_sigma / GWSAMPLESIZE;
    m_gs_inv = 1.0 / gw_stepsize;

    REAL s0 = 0;
    for (INT i = 0; i < GWSAMPLESIZE; i ++) {
        m_gw_lut.push_back(exp(-((gw_stepsize * i) *
                                 (gw_stepsize * i)) * gw_c0) * gw_c1);
        s0 += m_gw_lut.back();
    }
#else
    REAL c0 = m_tw_spacing * 0.5;
	REAL c1 = m_tw_sigma*sqrt(2.0);
	REAL gw_stepsize = TRUNCATION * m_tw_sigma / GWSAMPLESIZE;
	m_gs_inv = 1.0 / gw_stepsize;

	REAL s0 = 0;
	for (int i = 0; i < GWSAMPLESIZE; i ++) {
		m_gw_lut.push_back((erf((i*gw_stepsize + c0)/c1) - erf((i*gw_stepsize - c0)/c1)) / 2.0);
	    s0 += m_gw_lut.back();
	}
#endif    
    
#if 0    
    mexPrintf("Gaussian window: %d, %.10f, %.10f\n", GWSAMPLESIZE, s0, s0 / GWSAMPLESIZE);
#endif
#endif
    
}

ImageRayTracer::~ImageRayTracer()
{
}
    
inline bool ImageRayTracer::hitCheck(const REAL p0x, const REAL p0y, const REAL p0z,
                                     const REAL dx, const REAL dy, const REAL dz,
                                     REAL& t_min, REAL& t_max)
{
#if CFOV_ENABLED
    
//    REAL fov_radius = 58 * m_vox_size_i * 0.5;
    REAL fov_radius = (m_vdim_i - 1) * m_vox_size_i * 0.5;
    REAL dx2dy2 = dx * dx + dy * dy;
    REAL det = dx2dy2 * fov_radius * fov_radius -
                 (dx * p0y - dy * p0x) * (dx * p0y - dy * p0x);

    if (det > 0.0) {
        REAL sqr_det = sqrtf(det);
        REAL xdx_ydy = -(dx * p0x + dy * p0y);
        t_min = (xdx_ydy - sqr_det) / dx2dy2;
        t_max = (xdx_ydy + sqr_det) / dx2dy2;
        return true;
    } else {
        return false;
    }

#else    
    // these are tricks to avoid the `divide by zero' error
    REAL ddx = (dx == 0) ? 1e-20 : dx;
    REAL ddy = (dy == 0) ? 1e-20 : dy;
    REAL ddz = (dz == 0) ? 1e-20 : dz;

    // parameter (need change if ray goes from p1 to p0)
    REAL tx0 = (ddx > 0) ? ((m_vbd_x0 - p0x) / dx) : ((m_vbd_x1 - p0x) / ddx);
    REAL tx1 = (ddx > 0) ? ((m_vbd_x1 - p0x) / dx) : ((m_vbd_x0 - p0x) / ddx);
    REAL ty0 = (ddy > 0) ? ((m_vbd_y0 - p0y) / dy) : ((m_vbd_y1 - p0y) / ddy);
    REAL ty1 = (ddy > 0) ? ((m_vbd_y1 - p0y) / dy) : ((m_vbd_y0 - p0y) / ddy);
    REAL tz0 = (ddz > 0) ? ((m_vbd_z0 - p0z) / dz) : ((m_vbd_z1 - p0z) / ddz);
    REAL tz1 = (ddz > 0) ? ((m_vbd_z1 - p0z) / dz) : ((m_vbd_z0 - p0z) / ddz);
    /*
    if (dx == 0) {
        tx0 = -99.9;
        tx1 = +99.9;
    }

    if (dy == 0) {
        ty0 = -99.9;
        ty1 = +99.9;
    }

    if (dz == 0) {
        tz0 = -99.9;
        tz1 = +99.9;
    }     
     */
    
    // determine min and max
    t_min = std::max(std::max(tx0, std::max(ty0, tz0)), 0.0); 
    t_max = std::min(std::min(tx1, std::min(ty1, tz1)), 1.0);
    return (t_min < t_max);
#endif
}

inline void ImageRayTracer::calcTOFLOREndPoints(const INT tbin_id,
                                                REAL& p0x, REAL& p0y, REAL& p0z,
                                                REAL& p1x, REAL& p1y, REAL& p1z,
                                                REAL& nrm_dx, REAL& nrm_dy, REAL& nrm_dz)
{
    // ray direction vector
    REAL dx = p1x - p0x;
    REAL dy = p1y - p0y;
    REAL dz = p1z - p0z;

    // lor center
    REAL pcx = (p0x + p1x) * 0.5;
    REAL pcy = (p0y + p1y) * 0.5;
    REAL pcz = (p0z + p1z) * 0.5;

    // (inverse) length
    REAL inv_ray_length = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);

    // normalized dir
    nrm_dx = dx * inv_ray_length;
    nrm_dy = dy * inv_ray_length;
    nrm_dz = dz * inv_ray_length;

    REAL tbin_offset1 = tbin_id * m_tw_spacing - TRUNCATION * m_tw_sigma;
    REAL tbin_offset2 = tbin_id * m_tw_spacing + TRUNCATION * m_tw_sigma;

    // get new end-points
    p0x = pcx + tbin_offset1 * nrm_dx;
    p0y = pcy + tbin_offset1 * nrm_dy;
    p0z = pcz + tbin_offset1 * nrm_dz;
    p1x = pcx + tbin_offset2 * nrm_dx;
    p1y = pcy + tbin_offset2 * nrm_dy;
    p1z = pcz + tbin_offset2 * nrm_dz;
}

REAL ImageRayTracer::fproj(REAL p0x, REAL p0y, REAL p0z,
                           REAL p1x, REAL p1y, REAL p1z,
                           const IMAGE_DATA_TYPE* img, 
                           const INT tbin_id)
{
#ifdef USE_TOF
    // normalized dir
    REAL nrm_dx;
    REAL nrm_dy;
    REAL nrm_dz;
    calcTOFLOREndPoints(tbin_id, 
                        p0x, p0y, p0z, 
                        p1x, p1y, p1z, 
                        nrm_dx, nrm_dy, nrm_dz);
    // timing bin center coord
    REAL tbc_x = (p0x + p1x) * 0.5;
    REAL tbc_y = (p0y + p1y) * 0.5;
    REAL tbc_z = (p0z + p1z) * 0.5;
#endif
    
    REAL dx = p1x - p0x;
    REAL dy = p1y - p0y;
    REAL dz = p1z - p0z;
    REAL t_min, t_max;
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return 0;
    }
#ifdef USE_TOF

    t_min = std::max(t_min, 0.0);
    t_max = std::min(t_max, 1.0);

    // double check if the ray segment is still inside the FOV
    if (t_min > t_max) {
        return 0.0;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    REAL coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif    
    
    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz);     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the max index value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);
    
    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;
    
    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;

    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    REAL out = 0.0;
    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0;
        
#ifdef USE_TOF // TOF

#ifndef DIST_NO_SQRT
        REAL tdx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j - tbc_x; //m_x_cuts[j0] - tbc_x;
        REAL tdy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i - tbc_y; //m_y_cuts[m_img_size_i - i0 - 1] - tbc_y;
        REAL tdz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k - tbc_z; //m_z_cuts[k0] - tbc_z;
        REAL vox_to_tbc_dist = sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        REAL vx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j; //m_x_cuts[j0];
        REAL vy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i; //m_y_cuts[m_img_size_i - i0 - 1];
        REAL vz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k; //m_z_cuts[k0];
        REAL vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif

        INT gw_bin = INT(vox_to_tbc_dist * m_gs_inv);
        if (gw_bin < GWSAMPLESIZE) {
            out += img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] * ww * m_gw_lut[gw_bin];
        }

#else // nonTOF
        // put value to image domain
        // calc voxel index
        out += ww * img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0]; //img(i0, j0, k0);
#endif
        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

    return out;
}

void ImageRayTracer::bproj(REAL p0x, REAL p0y, REAL p0z,
                           REAL p1x, REAL p1y, REAL p1z,
                           const REAL weight, IMAGE_DATA_TYPE* img, 
                           const INT tbin_id)
{
    if (weight == 0) {
        return;
    }
#ifdef USE_TOF
    // normalized dir
    REAL nrm_dx;
    REAL nrm_dy;
    REAL nrm_dz;
    calcTOFLOREndPoints(tbin_id, p0x, p0y, p0z, p1x, p1y, p1z, nrm_dx, nrm_dy, nrm_dz);
    
    // timing bin center coord
    REAL tbc_x = (p0x + p1x) * 0.5f;
    REAL tbc_y = (p0y + p1y) * 0.5f;
    REAL tbc_z = (p0z + p1z) * 0.5f;
#endif
    
    REAL dx = p1x - p0x;
    REAL dy = p1y - p0y;
    REAL dz = p1z - p0z;
    REAL t_min, t_max;

    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return;
    }
    
#ifdef USE_TOF

    t_min = std::max(t_min, 0.0);
    t_max = std::min(t_max, 1.0);
    // double check if the ray segment is still inside the FOV
    if (t_min > t_max) {
        return;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    REAL coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif    
    
    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz) * weight;     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the index max value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);

    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;

    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;
    
    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
//        mexPrintf("%d %d %d, %.20f, tx=%.20f, ty=%.20f, tz=%.20f\n", i, j, k, tm0, tx,ty,tz);
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0; 

#ifdef USE_TOF // TOF
////<---
#ifndef DIST_NO_SQRT
        REAL tdx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j - tbc_x; //m_x_cuts[j0] - tbc_x;
        REAL tdy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i - tbc_y; //m_y_cuts[m_img_size_i - i0 - 1] - tbc_y;
        REAL tdz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k - tbc_z; //m_z_cuts[k0] - tbc_z;
        REAL vox_to_tbc_dist = sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        REAL vx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j; //m_x_cuts[j0];
        REAL vy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i; //m_y_cuts[m_img_size_i - i0 - 1];
        REAL vz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k; //m_z_cuts[k0];
        REAL vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif

#if 1        
        INT gw_bin = INT(vox_to_tbc_dist * m_gs_inv); 
        if (gw_bin < GWSAMPLESIZE) {
#ifdef USE_OMP
            #pragma omp atomic
#endif        
            img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] += ww * m_gw_lut[gw_bin];
        }
#endif

#else // nonTOF
        
        // put value to image domain
        // calc voxel index
#ifdef USE_OMP
//#if BP_ATOMIC
        #pragma omp atomic
//#endif
#endif        
        img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] += ww; 
#endif
        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

}

// A->S->B' v1
REAL ImageRayTracer::fproj_SSS_AS(REAL p0x, REAL p0y, REAL p0z, // A
                                  REAL sx, REAL sy, REAL sz, // S
                                  REAL virtual_x, REAL virtual_y, REAL virtual_z, // B'
                                  const IMAGE_DATA_TYPE* img, 
                                  const INT tbin_id)
{
#ifdef USE_TOF
    // normalized dir
    REAL nrm_dx;
    REAL nrm_dy;
    REAL nrm_dz;
    REAL ax0 = p0x;
    REAL ay0 = p0y;
    REAL az0 = p0z;
    REAL ax1 = virtual_x;
    REAL ay1 = virtual_y;
    REAL az1 = virtual_z;
    calcTOFLOREndPoints(tbin_id, 
                        ax0, ay0, az0, 
                        ax1, ay1, az1, 
                        nrm_dx, nrm_dy, nrm_dz);
    // timing bin center coord
    REAL tbc_x = (ax0 + ax1) * 0.5;
    REAL tbc_y = (ay0 + ay1) * 0.5;
    REAL tbc_z = (az0 + az1) * 0.5;
#if 0
mexPrintf("\t[%f,%f,%f]\n", tbc_x, tbc_y, tbc_z);
#endif

#endif
    
    REAL dx = sx - p0x;
    REAL dy = sy - p0y;
    REAL dz = sz - p0z;
    REAL t_min, t_max;
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return 0;
    }

#if 0
mexPrintf("\t[tmin:%f, tmax:%f]\n", t_min, t_max);
#endif

#ifdef USE_TOF

    t_min = std::max(t_min, 0.0);
    t_max = std::min(t_max, 1.0);

    // double check if the ray segment is still inside the FOV
    if (t_min > t_max) {
        return 0.0;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    REAL coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif    
    
    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz);     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the max index value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);
    
    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;
    
    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;

    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    REAL out = 0.0;
    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0;
        
#ifdef USE_TOF // TOF

#ifndef DIST_NO_SQRT
        REAL tdx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j - tbc_x; //m_x_cuts[j0] - tbc_x;
        REAL tdy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i - tbc_y; //m_y_cuts[m_img_size_i - i0 - 1] - tbc_y;
        REAL tdz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k - tbc_z; //m_z_cuts[k0] - tbc_z;
        REAL vox_to_tbc_dist = sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        REAL vx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j; //m_x_cuts[j0];
        REAL vy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i; //m_y_cuts[m_img_size_i - i0 - 1];
        REAL vz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k; //m_z_cuts[k0];
        REAL vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif

        INT gw_bin = INT(vox_to_tbc_dist * m_gs_inv);
        if (gw_bin < GWSAMPLESIZE) {
            out += img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] * ww * m_gw_lut[gw_bin];
        }

#else // nonTOF
        // put value to image domain
        // calc voxel index
        out += ww * img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0]; //img(i0, j0, k0);
#endif
        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

    return out;
}

void ImageRayTracer::calcTOFBinCenter(REAL p0x, REAL p0y, REAL p0z,
                                      REAL p1x, REAL p1y, REAL p1z,
                                      REAL& nrm_dx, REAL& nrm_dy, REAL& nrm_dz,
                                      const std::vector<INT>& tbin,
                                      std::vector<REAL>& tbin_center_list_x,
                                      std::vector<REAL>& tbin_center_list_y,
                                      std::vector<REAL>& tbin_center_list_z)
{
    // ray direction vector
    REAL dx = p1x - p0x;
    REAL dy = p1y - p0y;
    REAL dz = p1z - p0z;

    // lor center
    REAL pcx = (p0x + p1x) * 0.5;
    REAL pcy = (p0y + p1y) * 0.5;
    REAL pcz = (p0z + p1z) * 0.5;

    // (inverse) length
    REAL inv_ray_length = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);

    // normalized dir
    nrm_dx = dx * inv_ray_length;
    nrm_dy = dy * inv_ray_length;
    nrm_dz = dz * inv_ray_length;

    INT c = 0;
    for (INT nn=0; nn<tbin.size(); nn++) {

        INT n = tbin[nn];
        tbin_center_list_x[c] = pcx + n * m_tw_spacing * nrm_dx;
        tbin_center_list_y[c] = pcy + n * m_tw_spacing * nrm_dy;
        tbin_center_list_z[c] = pcz + n * m_tw_spacing * nrm_dz;
        c ++;
    }

}

// A->S->B' v2
void ImageRayTracer::fproj_TOF_SSS_AS(REAL p0x, REAL p0y, REAL p0z, // A
                                      REAL sx, REAL sy, REAL sz, // S
                                      REAL virtual_x, REAL virtual_y, REAL virtual_z, // B'
                                      const IMAGE_DATA_TYPE* img, 
                                      const std::vector<INT>& tbin,
                                      REAL* out)
{
    INT nbin = tbin.size(); //tbin_max * 2 + 1;
    // initial values
    memset(out, 0, sizeof(REAL) * nbin);

    REAL nrm_dx;
    REAL nrm_dy;
    REAL nrm_dz;
    REAL ax0 = p0x;
    REAL ay0 = p0y;
    REAL az0 = p0z;
    REAL ax1 = virtual_x;
    REAL ay1 = virtual_y;
    REAL az1 = virtual_z;
    std::vector<REAL> tbin_center_list_x(nbin);
    std::vector<REAL> tbin_center_list_y(nbin);
    std::vector<REAL> tbin_center_list_z(nbin);
    calcTOFBinCenter(ax0, ay0, az0, ax1, ay1, az1,
                     nrm_dx, nrm_dy, nrm_dz, tbin,
                     tbin_center_list_x,
                     tbin_center_list_y,
                     tbin_center_list_z);

    REAL dx = sx - p0x;
    REAL dy = sy - p0y;
    REAL dz = sz - p0z;
    REAL t_min, t_max;
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return;
    }

    // this is TOF
    t_min = std::max(t_min, 0.0);
    t_max = std::min(t_max, 1.0);

    // double check if the ray segment is still inside the FOV
    if (t_min > t_max) {
        return;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    std::vector<REAL> coef_d_list(nbin);
    for (int i = 0; i < nbin; i ++) {
        coef_d_list[i] = -(nrm_dx * tbin_center_list_x[i] + 
                           nrm_dy * tbin_center_list_y[i] + 
                           nrm_dz * tbin_center_list_z[i]);
    }
#endif

    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz);     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the max index value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);
    
    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;
    
    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;

    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0;
        
#ifndef DIST_NO_SQRT
        for (INT tt = 0; tt < nbin; tt ++) {
            REAL tdx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j - tbin_center_list_x[tt]; 
            REAL tdy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i - tbin_center_list_y[tt]; 
            REAL tdz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k - tbin_center_list_z[tt]; 
            REAL vox_to_tbc_dist = sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        REAL vx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j; //m_x_cuts[j0];
        REAL vy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i; //m_y_cuts[m_img_size_i - i0 - 1];
        REAL vz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k; //m_z_cuts[k0];
        for (INT tt = 0; tt < nbin; tt ++) {
            REAL vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d_list[tt]);
#endif
            INT gw_bin = INT(vox_to_tbc_dist * m_gs_inv);
            if (gw_bin < GWSAMPLESIZE) {
                out[tt] += img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] * ww * m_gw_lut[gw_bin];
            }
        }

        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

}

// A'->S->B
REAL ImageRayTracer::fproj_SSS_SB(REAL virtual_x, REAL virtual_y, REAL virtual_z, // A
                                  REAL sx, REAL sy, REAL sz, // S
                                  REAL p1x, REAL p1y, REAL p1z, // B'
                                  const IMAGE_DATA_TYPE* img, 
                                  const INT tbin_id)
{
#ifdef USE_TOF
    // normalized dir
    REAL nrm_dx;
    REAL nrm_dy;
    REAL nrm_dz;
    REAL ax0 = virtual_x;
    REAL ay0 = virtual_y;
    REAL az0 = virtual_z;
    REAL ax1 = p1x;
    REAL ay1 = p1y;
    REAL az1 = p1z;
    calcTOFLOREndPoints(tbin_id, 
                        ax0, ay0, az0, 
                        ax1, ay1, az1, 
                        nrm_dx, nrm_dy, nrm_dz);
    // timing bin center coord
    REAL tbc_x = (ax0 + ax1) * 0.5;
    REAL tbc_y = (ay0 + ay1) * 0.5;
    REAL tbc_z = (az0 + az1) * 0.5;
#endif
    
    REAL p0x = sx;
    REAL p0y = sy;
    REAL p0z = sz;

    REAL dx = p1x - sx;
    REAL dy = p1y - sy;
    REAL dz = p1z - sz;

    REAL t_min, t_max;
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return 0;
    }
#ifdef USE_TOF

    t_min = std::max(t_min, 0.0);
    t_max = std::min(t_max, 1.0);

    // double check if the ray segment is still inside the FOV
    if (t_min > t_max) {
        return 0.0;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    REAL coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif    
    
    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz);     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the max index value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);
    
    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;
    
    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;

    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    REAL out = 0.0;
    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0;
        
#ifdef USE_TOF // TOF

#ifndef DIST_NO_SQRT
        REAL tdx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j - tbc_x; //m_x_cuts[j0] - tbc_x;
        REAL tdy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i - tbc_y; //m_y_cuts[m_img_size_i - i0 - 1] - tbc_y;
        REAL tdz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k - tbc_z; //m_z_cuts[k0] - tbc_z;
        REAL vox_to_tbc_dist = sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        REAL vx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j; //m_x_cuts[j0];
        REAL vy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i; //m_y_cuts[m_img_size_i - i0 - 1];
        REAL vz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k; //m_z_cuts[k0];
        REAL vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif

        INT gw_bin = INT(vox_to_tbc_dist * m_gs_inv);
        if (gw_bin < GWSAMPLESIZE) {
            out += img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] * ww * m_gw_lut[gw_bin];
        }

#else // nonTOF
        // put value to image domain
        // calc voxel index
        out += ww * img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0]; //img(i0, j0, k0);
#endif
        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

    return out;
}

// A'->S->B v2
void ImageRayTracer::fproj_TOF_SSS_SB(REAL virtual_x, REAL virtual_y, REAL virtual_z, // A'
                                      REAL sx, REAL sy, REAL sz, // S
                                      REAL p1x, REAL p1y, REAL p1z, // B
                                      const IMAGE_DATA_TYPE* img, 
                                      const std::vector<INT>& tbin,
                                      REAL* out)
{
    INT nbin = tbin.size(); //tbin_max * 2 + 1;
    // initial values
    memset(out, 0, sizeof(REAL) * nbin);

    REAL nrm_dx;
    REAL nrm_dy;
    REAL nrm_dz;
    REAL ax0 = virtual_x;
    REAL ay0 = virtual_y;
    REAL az0 = virtual_z;
    REAL ax1 = p1x;
    REAL ay1 = p1y;
    REAL az1 = p1z;
    std::vector<REAL> tbin_center_list_x(nbin);
    std::vector<REAL> tbin_center_list_y(nbin);
    std::vector<REAL> tbin_center_list_z(nbin);
    calcTOFBinCenter(ax0, ay0, az0, ax1, ay1, az1,
                     nrm_dx, nrm_dy, nrm_dz, tbin,
                     tbin_center_list_x,
                     tbin_center_list_y,
                     tbin_center_list_z);
    
    REAL p0x = sx;
    REAL p0y = sy;
    REAL p0z = sz;
    REAL dx = p1x - sx;
    REAL dy = p1y - sy;
    REAL dz = p1z - sz;

    REAL t_min, t_max;
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return;
    }

    t_min = std::max(t_min, 0.0);
    t_max = std::min(t_max, 1.0);

    // double check if the ray segment is still inside the FOV
    if (t_min > t_max) {
        return;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    std::vector<REAL> coef_d_list(nbin);
    for (int i = 0; i < nbin; i ++) {
        coef_d_list[i] = -(nrm_dx * tbin_center_list_x[i] + 
                           nrm_dy * tbin_center_list_y[i] + 
                           nrm_dz * tbin_center_list_z[i]);
    }
#endif

    
    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz);     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the max index value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);
    
    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;
    
    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;

    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0;
        
#ifndef DIST_NO_SQRT
        for (INT tt = 0; tt < nbin; tt ++) {
            REAL tdx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j - tbin_center_list_x[tt]; 
            REAL tdy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i - tbin_center_list_y[tt]; 
            REAL tdz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k - tbin_center_list_z[tt]; 
            REAL vox_to_tbc_dist = sqrt(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        REAL vx = m_vbd_x0 + (j0 + 0.5) * m_vox_size_j; //m_x_cuts[j0];
        REAL vy = m_vbd_y1 - (i0 + 0.5) * m_vox_size_i; //m_y_cuts[m_img_size_i - i0 - 1];
        REAL vz = m_vbd_z0 + (k0 + 0.5) * m_vox_size_k; //m_z_cuts[k0];
        for (INT tt = 0; tt < nbin; tt ++) {
            REAL vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d_list[tt]);
#endif

            INT gw_bin = INT(vox_to_tbc_dist * m_gs_inv);
            if (gw_bin < GWSAMPLESIZE) {
                out[tt] += img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0] * ww * m_gw_lut[gw_bin];
            }
        }

        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

}

REAL ImageRayTracer::fproj_nonTOF(REAL p0x, REAL p0y, REAL p0z,
                                  REAL p1x, REAL p1y, REAL p1z,
                                  const IMAGE_DATA_TYPE* img)
{    
    REAL dx = p1x - p0x;
    REAL dy = p1y - p0y;
    REAL dz = p1z - p0z;
    REAL t_min, t_max;
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, t_min, t_max)) {
        return 0;
    }
    
    // calc length of ray
    REAL w0 = sqrt(dx * dx + dy * dy + dz * dz);     
    REAL p1stx = p0x + dx * t_min;
    REAL p1sty = p0y + dy * t_min;
    REAL p1stz = p0z + dz * t_min;
    
    // see the definition of coordination system
    // index of the first voxel hit by ray
    INT j = INT((p1stx - m_vbd_x0) / m_vox_size_j); // 1: x
    j = std::min(j, m_vdim_j - 1); // because the max index value is m_vdim_*-1

    INT i = INT((m_vbd_y1 - p1sty) / m_vox_size_i); // 0: y
    i = std::min(i, m_vdim_i - 1);

    INT k = INT((p1stz - m_vbd_z0) / m_vox_size_k); // 2 : z
    k = std::min(k, m_vdim_k - 1);
    
    // initial boundary
    REAL bx0 = (dx > 0) ? m_vbd_x0 + j * m_vox_size_j : m_vbd_x0 + (j + 1) * m_vox_size_j;
    REAL by0 = (dy > 0) ? m_vbd_y0 + (m_vdim_i - i - 1) * m_vox_size_i : 
                          m_vbd_y0 + (m_vdim_i - i) * m_vox_size_i;
    REAL bz0 = (dz > 0) ? m_vbd_z0 + k * m_vox_size_k : m_vbd_z0 + (k + 1) * m_vox_size_k;

    // step for update index
    INT di = (dy > 0) ? 1 : -1;
    INT dj = (dx > 0) ? 1 : -1;
    INT dk = (dz > 0) ? 1 : -1;
    
    // step for forward ray
    REAL ddx = (dx > 0) ? m_vox_size_j : -m_vox_size_j;
    REAL ddy = (dy > 0) ? m_vox_size_i : -m_vox_size_i;
    REAL ddz = (dz > 0) ? m_vox_size_k : -m_vox_size_k;

    INT idx;
    INT m = 0; // counter for nonzero element number
    REAL tm0 = t_min, tm1;
    REAL tx = (fabs(dx) < FLTMIN) ? 99.9 : (bx0 + ddx - p0x) / dx;
    REAL ty = (fabs(dy) < FLTMIN) ? 99.9 : (by0 + ddy - p0y) / dy;
    REAL tz = (fabs(dz) < FLTMIN) ? 99.9 : (bz0 + ddz - p0z) / dz;
    REAL dtx = (fabs(dx) < FLTMIN) ? 99.9 : ddx / dx;
    REAL dty = (fabs(dy) < FLTMIN) ? 99.9 : ddy / dy;
    REAL dtz = (fabs(dz) < FLTMIN) ? 99.9 : ddz / dz;

    REAL out = 0.0;
    // start main loop
    do {

        // current voxel subscripts
        INT i0 = i;
        INT j0 = j;
        INT k0 = k;
        
        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += dtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += dty;
            } else {
                tm1 = tz;
                k += dk;
                tz += dtz;
            }
        }

        if (tm1 > t_max) {
            tm1 = t_max;
        }

        // calc weight
        REAL ww = (tm1 - tm0) * w0;        
        // put value to image domain
        // calc voxel index
        out += ww * img[k0 * m_vdim_ixj + j0 * m_vdim_i + i0]; //img(i0, j0, k0);
        tm0 = tm1; // save previous result
    } while (fabs(tm0 - t_max) > FLTMIN);

    return out;
}

#endif
