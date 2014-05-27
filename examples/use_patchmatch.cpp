/*
 #
 #  File        : patchmatch.h
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

#include <iostream>
#define cimg_plugin "plugins/patchmatch.h"
#include "CImg.h"

using namespace cimg_library;
using namespace std;

int main(int argc, char **argv) {
  cimg_usage("PatchMatch: approximate nearest neighbour patch matching between two images img0 and img1.");
  cimg_help("Input / Output parameters:\n--------------------------");
  const char *file0 = cimg_option("-i0", (char*)0, "Image 0.");
  const char *file1 = cimg_option("-i1", (char*)0, "Image 1.");
  const int P = cimg_option("-s", 15, "Patch size.");
  const int N = cimg_option("-n", 5, "Number of iterations.");

  if (!file0 || !file1) {
    cerr<<"Must specify images with -i0 and -i1 options (-h for help).\n";
    exit(0);
  }
  srand(time(0));

  CImg<int> img0(file0);
  CImg<int> img1(file1);

  CImg<int> off;
  off.patchMatch(img0, img1, P, N);

  off.get_vizFlow(50).display();
  CImg<int> imgrec(img0.get_reconstruct(img1, off));
  (img0, imgrec, img1).display();
}
