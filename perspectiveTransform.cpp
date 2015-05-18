/*
 * perspectiveTransform (c) 2008-2015 Luc Deschenaux
 *
 * Thanks to pseudocode from developpez.net for the forward mapping
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Additional Terms:
 *
 *      You are required to preserve legal notices and author attributions in
 *      that material or in the Appropriate Legal Notices displayed by works
 *      containing it.
 *
 */

#define __perspectiveTransform
#include "perspectiveTransform.hpp"

/**
 * perspectiveTransformImage()
 *
 * Specified region in source image is warped to fit destination image
 * When map is not null, it is a pointer to an array of dimension dst->width()*dst->height()*siezof(ptPoint)
 * that will contain the original pixel positions
 *
 * Points must be specified in clockwise order (top-left, top,right, bot-right, bot-left)
 *
 * @param src   source image
 * @param P     point list
 * @param dst   target image
 * @param map   optional pixel map
 * @param H     homography matrix
*/

void perspectiveTransformImage(CImg<unsigned char> *src,CImg<unsigned char> *dst,ptPoint *map,double *H) {

  double x;
  double y;
  double width=(double)dst->width();
  double height=(double)dst->height();
  double *S;

  // pour chaque pixel de l'image cible  
  for(y=0;y<height;y++) {
    for(x=0;x<width;x++) {
 
      // conversion dans le repère orthonormé (u,v) [0,1]x[0,1]
      double u = x/width;
      double v = y/height;

      // passage dans le repère perspective
      projective_mapping(&u, &v, H);

      // coordonées dans l'image source
      double sx=u;
      double sy=v;
    
      // read source pixel value 
      unsigned char r=src->cubic_atXY(sx,sy,0,0);
      unsigned char g=src->cubic_atXY(sx,sy,0,1);
      unsigned char b=src->cubic_atXY(sx,sy,0,2);
      unsigned char a=src->cubic_atXY(sx,sy,0,3);
    
      // write destination pixel value  
      int dx=x;
      int dy=y;
      dst[0](dx,dy,0,0)=r;
      dst[0](dx,dy,0,1)=g;
      dst[0](dx,dy,0,2)=b;
      dst[0](dx,dy,0,3)=a;

      // store original pixel position
      if (map) {
        size_t index=dy*width+dx;
        map[index].x=sx;
        map[index].y=sy;
      }
    }
  }
}

/**
 * inversePerspectiveTransformImage()
 *
 * Source image is warped to fit the specified region in dest image
 * When map is not null, it is a pointer to an array of dimension dst->width()*dst->height()*siezof(ptPoint)
 * that will contain the original pixel positions
 *
 * @param src                     source image
 * @param dst                     target image
 * @param map                     optional pixel map
 * @param H                       homography matrix
 * @param transparentPixelsOnly   if set, this flag restrict 
*/

void inversePerspectiveTransformImage(CImg<unsigned char> *src,CImg<unsigned char> *dst,ptPoint *map,double *H, bool transparentPixelsOnly) {

  double width=dst->width();
  double height=dst->height();
  double *S;

  // for every pixel in the target image
  for(double dy=0;dy<height;dy+=1) {
    for(double dx=0;dx<width;dx+=1) {

      // skip non transparent pixels if requested
      if (transparentPixelsOnly && dst[0](dx,dy,0,3)>0) {
        continue;
      }

      // coordonnées du pixel de destination
      double u = dx;
      double v = dy;

      // passage dans le repère orthonormé
      inverse_projective_mapping(&u, &v, H+9);

      // coordonnées dans l'image source
      double sx=u*src->width();
      double sy=v*src->height();

      // read source pixel value 
      unsigned char r=src->cubic_atXY(sx,sy,0,0);
      unsigned char g=src->cubic_atXY(sx,sy,0,1);
      unsigned char b=src->cubic_atXY(sx,sy,0,2);
      unsigned char a=src->cubic_atXY(sx,sy,0,3);
    
      // write destination pixel value  
      dst[0](dx,dy,0,0)=r;
      dst[0](dx,dy,0,1)=g;
      dst[0](dx,dy,0,2)=b;
      dst[0](dx,dy,0,3)=a;

      // store original pixel position
      if (map) {
        size_t index=dy*width+dx;
        map[index].x=sx;
        map[index].y=sy;
      }
    }
  }

} // inversePerspectiveTransformImage

void perspectiveTransform(CImg<unsigned char> *dst, ptPoint *I, double *H) {

  // conversion dans le repère orthonormé (u,v) [0,1]x[0,1]
  double u = (double)I->x/(double)dst->width();
  double v = (double)I->y/(double)dst->height();

  // passage dans le repère perspective
  projective_mapping(&u, &v, H);

  // pixel correspondant dans le quadrilatère
  I->x=u;
  I->y=v;

} // perspectiveTransform

void inversePerspectiveTransform(CImg<unsigned char> *src, ptPoint *I, double *H) {

  double u=I->x;
  double v=I->y;

  // passage dans le repère orthonormé
  inverse_projective_mapping(&u, &v, H+9);

  // pixel correspondant dans le quadrilatère
  I->x=u*(double)src->width();
  I->y=v*(double)src->height();

} // inversePerspectiveTransform

/*
 The algorithm is coming from:
 Projective Mappings for Image Warping (1989)
 by Paul Heckbert
 in Fundamentals of Texture Mapping and Image Warping (Paul Heckbert,
 Master’s Thesis), U.C.Berkeley

 http://www.cs.utah.edu/classes/cs6964-whitaker/heckbert_projective.pdf 

 The case of parallelograms is not handled.

 You can combine the transformations with matricial product to map one
 quadrilateral to another one via the square.

 Thanks to pseudocode from developpez.net for the forward mapping

 Points must be specified in clockwise order (top-left, top,right, bot-right, bot-left)

*/

double *getPerspectiveTransform(ptPoint *P) {

  double *H=(double*)calloc(18,sizeof(double));
  double *adj=H+9;

  double sx = (P[0].x-P[1].x)+(P[2].x-P[3].x);
  double sy = (P[0].y-P[1].y)+(P[2].y-P[3].y);
  double dx1 = P[1].x-P[2].x;
  double dx2 = P[3].x-P[2].x;
  double dy1 = P[1].y-P[2].y;
  double dy2 = P[3].y-P[2].y;
 
  double z = (dx1*dy2)-(dy1*dx2);
  double g = ((sx*dy2)-(sy*dx2))/z;
  double h = ((sy*dx1)-(sx*dy1))/z;
 
  // matrice de transformation
  double a=H[0]=P[1].x-P[0].x+g*P[1].x;
  double b=H[1]=P[3].x-P[0].x+h*P[3].x;
  double c=H[2]=P[0].x;
  double d=H[3]=P[1].y-P[0].y+g*P[1].y;
  double e=H[4]=P[3].y-P[0].y+h*P[3].y;
  double f=H[5]=P[0].y;
  H[6]=g;
  H[7]=h;
  H[8]=1;

  // calcul de transformation inverse (matrice adjointe)
  adj[0]=e-f*h;
  adj[1]=c*h-b;
  adj[2]=b*f-c*e;
  adj[3]=f*g-d;
  adj[4]=a-c*g;
  adj[5]=c*d-a*f;
  adj[6]=d*h-e*g;
  adj[7]=b*g-a*h;
  adj[8]=a*e-b*d;

//  printf("a %f b %f c %f d %f e %f f %f g %f h %f i 1\n",a,b,c,d,e,f,g,h);
//  printf("adj: a %f b %f c %f d %f e %f f %f g %f h %f i 1\n",adj[0],adj[1],adj[2],adj[3],adj[4],adj[5],adj[6],adj[7]);
  return H;
}

void projective_mapping(double *u, double *v, double *H) {
  double x = (H[0]*u[0]+H[1]*v[0]+H[2])/(H[6]*u[0]+H[7]*v[0]+1); 
  double y = (H[3]*u[0]+H[4]*v[0]+H[5])/(H[6]*u[0]+H[7]*v[0]+1);
  u[0]=x;
  v[0]=y;
}

void inverse_projective_mapping(double *u, double *v, double *H) {
  double x = (H[0]*u[0]+H[1]*v[0]+H[2])/(H[6]*u[0]+H[7]*v[0]+1); 
  double y = (H[3]*u[0]+H[4]*v[0]+H[5])/(H[6]*u[0]+H[7]*v[0]+1);
  double z = (H[6]*u[0]+H[7]*v[0]+H[8])/(H[6]*u[0]+H[7]*v[0]+1);
  u[0]=x/z;
  v[0]=y/z;
}


