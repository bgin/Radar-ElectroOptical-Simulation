

#ifndef __GMS_RCS_ELLIPS_OGIVE_ZMM16R4_HPP__
#define __GMS_RCS_ELLIPS_OGIVE_ZMM16R4_HPP__ 010320230919


/*MIT License
Copyright (c) 2020 Bernard Gingold
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

namespace file_version {

    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM16R4_MAJOR = 1U;
    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM16R4_MINOR = 0U;
    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM16R4_MICRO = 0U;
    const unsigned int GMS_RCS_ELLIPS_OGIVE_ZMM16R4_FULLVER =
      1000U*GMS_RCS_ELLIPS_OGIVE_ZMM16R4_MAJOR+
      100U*GMS_RCS_ELLIPS_OGIVE_ZMM16R4_MINOR+
      10U*GMS_RCS_ELLIPS_OGIVE_ZMM16R4_MICRO;
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM16R4_CREATION_DATE = "01-03-2023 09:18 AM +00200 (WED 01 MAR 2023 GMT+2)";
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM16R4_BUILD_DATE    = __DATE__ ":" __TIME__;
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM16R4_AUTHOR        = "Programmer: Bernard Gingold, contact: beniekg@gmail.com";
    const char * const GMS_RCS_ELLIPS_OGIVE_ZMM16R4_DESCRIPTION   = "AVX512 optimized Ellipsoid and Ogive Radar Cross Section (analytic) functionality.";

}

#include <cstdint>
#include <immintrin.h>
#include "GMS_config.h"
#include "GMS_sleefsimdsp.hpp"


namespace gms {


          namespace radiolocation {


                      /*
                            High-frequency cross-section of perfectly
                            conducting ellipsoid.
                            Formula 5.1.54
                        */


                  

       }



}














#endif /*__GMS_RCS_ELLIPS_OGIVE_ZMM16R4_HPP__*/
