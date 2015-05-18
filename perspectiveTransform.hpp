
#ifndef __perspectiveTransform__
#define __perspectiveTransform__

#include <stdlib.h>
#include <math.h>

#ifndef __perspectiveTransform
#define __extern extern
#else
#define __extern
#endif

#ifndef cimg_version
#define cimg_use_tiff
#define cimg_display 0
#include <CImg.h>
using namespace cimg_library;
#endif         

typedef struct ptPoint {
  float x;
  float y;
} ptPoint;

// public
__extern double *getPerspectiveTransform(ptPoint *P);
__extern void perspectiveTransformImage(CImg<unsigned char> *src,CImg<unsigned char> *dst,ptPoint *map,double *H);
__extern void perspectiveTransform(CImg<unsigned char> *img,ptPoint *I,double *H);
__extern void inversePerspectiveTransform(CImg<unsigned char> *img,ptPoint *I,double *H);
__extern void inversePerspectiveTransformImage(CImg<unsigned char> *src,CImg<unsigned char> *dst,ptPoint *map,double *H,bool transparentPixelsOnly);

// private
__extern void projective_mapping(double *u, double *v, double *H);
__extern void inverse_projective_mapping(double *x, double *y, double *H);

#undef __extern

#endif

