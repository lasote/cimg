/*
 #
 #  File        : PatchMatch_plugin.h
 #                ( C++ header file - CImg plug-in )
 #
 #  Description : Plugin implementing the Patch Match algorithm to use
 #                with the CImg library.
 #                ( http://cimg.sourceforge.net )
 #
 #  Copyright   : Olivier D'Hondt
 #                (https://sites.google.com/site/dhondtolivier/)
 #
 #  License     : CeCILL v2.0
 #                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed by the CeCILL  license under French law and
 #  abiding by the rules of distribution of free software.  You can  use,
 #  modify and/ or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
 #
*/

// Visualize optical flow maps with HSV color coding.
CImg<float> get_vizFlow(float cutVal = 0)
{
  CImg<float> res((*this), "x,y,1,3", 0.0);

  // Normalizing offset magnitude
  CImg<float> mag((*this), "xy", 0.0);
  cimg_forXY(*this, x, y){
    mag(x, y) = std::sqrt((*this)(x, y, 0, 0)*(*this)(x, y, 0, 0) + (*this)(x , y, 0, 1)*(*this)(x, y, 0, 1));
  }

  if(cutVal)
    mag.cut(0, cutVal);

  mag /= mag.max();

  // Filling HSV values
  cimg_forXY(*this, x, y){
    float xx = -(*this)(x, y, 0, 0);
    float yy = -(*this)(x, y, 0, 1);

    float H = cimg::max(180*((std::atan2(yy, xx)/M_PI)+1.0), 0.0);
    float S = mag(x, y);
    float V = 1.0;

    res(x, y, 0, 0) = H;
    res(x, y, 0, 1) = S;
    res(x, y, 0, 2) = V;

  }

  res.HSVtoRGB();

  return res;
}

// Distance between two patches
T distPatch(const CImg<T> &img0, const CImg<T> &img1,
    int x0, int y0,
    int x1, int y1,
    int pSize)
{
  T d2 = 0;
  for(int c = 0; c < img0.spectrum(); c++){
    for(int y = 0; y < pSize; y++){
      for(int x = 0; x < pSize; x++){
        T d = (img0(x0 + x, y0 + y, 0, c) - img1(x1 + x, y1 + y, 0, c));
        d2 += d*d;
      }
    }
  }
  return d2;
}

// Patch match algorithm
template<typename Tt>
CImg<T> & patchMatch(const CImg<Tt> &img0, const CImg<Tt> &img1,
                    int patchSize, int nIter = 2, CImgDisplay* disp=NULL)
{
  if(img0.spectrum() != img1.spectrum())
    throw CImgInstanceException("Images must have the same number of channels.");

  if(!patchSize % 2){
    patchSize++;
    cimg::warn("Input patch size is even, adding 1.");
  }
  int w0 = img0.width();
  int h0 = img0.height();
  int w1 = img1.width();
  int h1 = img1.height();
  int nChannels = img0.spectrum();

  CImg<Tt> imgrec; // used only for display purpose

  if(disp) imgrec.assign(img0, "xy1c", 0.0);

  int P = patchSize;
  int H = P/2;

  // Zero padding borders
  CImg<Tt> img0big(w0+2*H, h0+2*H, 1, nChannels, 0);
  CImg<Tt> img1big(w1+2*H, h1+2*H, 1, nChannels, 0);

  // Try to penalize border patches
  img0big.rand(0,255);
  img1big.rand(0,255);

  img0big.draw_image(H, H, 0, 0, img0);
  img1big.draw_image(H, H, 0, 0, img1);

  CImg<T> off(w0, h0, 1, 2, 0);
  CImg<Tt> minDist(w0, h0, 1, 1, 0);

  // Initialize with random offsets
  cimg_forXY(off, x0, y0){
    int x1 = ((w1-1) * cimg::rand());
    int y1 = ((h1-1) * cimg::rand());
    off(x0, y0, 0, 0) = x1-x0;
    off(x0, y0, 0, 1) = y1-y0;
    minDist(x0, y0) = distPatch(img0big, img1big, x0, y0, x1, y1, P);
  }

  int xStart, yStart, xFinish, yFinish;
  int inc;
  for(int n = 0; n < nIter; n++){

    std::fprintf(stderr,"Iteration %d\n",n+1);
    // at odd iterations, reverse scan order
    if(n%2 == 0){
      xStart = 1; yStart = 1;
      xFinish = w0;
      yFinish = h0;
      inc = 1;
    }else{
      xStart = w0-2;
      yStart = h0-2;
      xFinish = -1; yFinish = -1;
      inc = -1;
    }
    for(int y = yStart; y != yFinish; y=y+inc)
      for(int x = xStart; x != xFinish; x=x+inc){
        // Propagation
        Tt d2 = 0.0;
        int x1 = x+off(x-inc, y, 0, 0);
        int y1 = y+off(x-inc, y, 0, 1);
        if(x1 >= 0 && x1 < w1 && y1 >= 0 && y1 < h1){ // propagate only if inside img1 bounds
          d2 = distPatch(img0big, img1big, x, y, x1, y1, P);
          if(d2<minDist(x, y)){
            minDist(x, y) = d2;
            off(x, y, 0, 0) = off(x-inc, y, 0, 0);
            off(x, y, 0, 1) = off(x-inc, y, 0, 1);
          }
        }
        x1 = x+off(x, y-inc, 0, 0);
        y1 = y+off(x, y-inc, 0, 1);
        if(x1 >= 0 && x1 < w1 && y1 >= 0 && y1 < h1){ // propagate only if inside img1 bounds
          d2 = distPatch(img0big, img1big, x, y, x1, y1, P);
          if(d2<minDist(x, y)){
            minDist(x, y) = d2;
            off(x, y, 0, 0) = off(x, y-inc, 0, 0);
            off(x, y, 0, 1) = off(x, y-inc, 0, 1);
          }
        }

        // Randomized search
        int wSizX = w1-1;
        int wSizY = h1-1;
        T offXCurr = off(x, y, 0, 0);
        T offYCurr = off(x, y, 0, 1);
        do{
          int wMinX = cimg::max(0, x+offXCurr-wSizX/2);
          int wMaxX = cimg::min(w1-1, x+offXCurr+wSizX/2);
          x1 = (wMaxX-wMinX) * cimg::rand() + wMinX;

          int wMinY = cimg::max(0, y+offYCurr-wSizY/2);
          int wMaxY = cimg::min(h1-1, y+offYCurr+wSizY/2);
          y1 = (wMaxY-wMinY) * cimg::rand() + wMinY;

          d2 = distPatch(img0big, img1big, x, y, x1, y1, P);

          if(d2 < minDist(x, y)){
            minDist(x, y) = d2;
            off(x, y, 0, 0) = x1-x;
            off(x, y, 0, 1) = y1-y;
          }
          wSizX /= 2;
          wSizY /= 2;
        }while(wSizX >= 1 && wSizY >= 1);
        // If a pointer to a CImgDisplay is passed as the last argument
        // the output of the algorithm is displayed as an animation
        // !! It slows down the algorithm a lot !!
        if(disp){
          if(x%(w0-1)==0){
            disp->display(
            (img0,
             imgrec.reconstruct(img1, off).get_normalize(0,255),
             off.get_vizFlow(100),
             img1));
          }
        }


      }
  }
  return off.move_to(*this);
}

// Reconstruct an image from an offset map and a query image
template<typename Tt>
CImg<T> & reconstruct(const CImg<T> &qimg, const CImg<Tt> &off)
{
  if((*this).spectrum() != qimg.spectrum())
    throw CImgInstanceException("Images must have the same number of channels.");
  if((*this).width() != off.width() || (*this).height() != off.height())
    throw CImgInstanceException("Offset map must have the same dimensions as input image.");
  cimg_forXY(off, x, y){
    int qx = x + off(x, y, 0, 0);
    int qy = y + off(x, y, 0, 1);
    cimg_forC(qimg, c){
      (*this)(x, y, 0, c) = qimg(qx, qy, 0, c);
    }
  }
  return (*this);
}

template<typename Tt>
CImg<T> get_reconstruct(const CImg<T> &qimg, const CImg<Tt> &off)
{
  return CImg<T>(*this, false).reconstruct(qimg, off);
}
