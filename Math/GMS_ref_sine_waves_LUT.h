
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

#ifndef __GMS_REF_SINE_WAVES_LUT_H__
#define __GMS_REF_SINE_WAVES_LUT_H__ 120920250813



/* Reference pure sine waves LUT arrays.
   Valid for Amplitude=1, Period=1, Samples per Period (32-40)
*/

#if !defined(REF_SINE_WAVES_LUT_USE_REP_HEX)
#define REF_SINE_WAVES_LUT_USE_REP_HEX 0
#endif 

namespace gms
{

namespace math 
{

// sin(2Pi*t/T) [amp=1,period=1,spp=32], spp i.e. samples per period
extern const float LUT_sin_1_1_32[32];

// sin(2Pi*t/T) [amp=1,period=1,samp=64]
extern const float LUT_sin_1_1_64[64];

// sin(2Pi*t/T) [amp=1,period=1,samp=128]
extern const float LUT_sin_1_1_128[128];

// sin(2Pi*t/T) [amp=1,period=1,samp=256]
extern const float LUT_sin_1_1_256[256];

// sin(2Pi*t/T) [amp=1,period=1,samp=512]
extern const float LUT_sin_1_1_512[512];

// sin(2Pi*t/T) [amp=1,period=1,samp=1024]
extern const float LUT_sin_1_1_1024[1024];

// sin(2Pi*t/T) [amp=1,period=1,samp=2048]
extern const float LUT_sin_1_1_2048[2048];

// sin(2Pi*t/T) [amp=1,period=1,samp=4096]
extern const float LUT_sin_1_1_4096[4096];


} // math 

} // gms

















#endif /*__GMS_REF_SINE_WAVES_LUT*/