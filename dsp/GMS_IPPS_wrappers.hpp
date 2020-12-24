
#ifndef __GMS_IPPS_WRAPPERS_HPP__
#define __GMS_IPPS_WRAPPERS_HPP__

#include <cstdint>
#include <stdio.h>
#include <ipps.h>
#include <ipp.h>
#include <ippcore.h>
#include "GMS_config.h"



__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
static inline
void gms_ippsGetLibVersion() {

     const IppLibraryVersion * __restrict ptrLib = NULL;
     ptrLib = ippsGetLibVersion();
     printf("major = %d\n",lib->major);
     printf("minor = %d\n",lib->minor);
     printf("majorBuild = %d\n",lib->majorBuild);
     printf("build = %d\n",lib->build);
     printf("targetCpu = %c%c%c%c\n",lib->targetCpu[0],lib->targetCpu[1],lib->targetCpu[2],
             lib->targetCpu[3]);
     printf("Name = %s\n", lib->Name);
     printf("Version = %s\n", lib->Version);
     printf("BuildDate = %s\n", lib->BuildDate);
}

__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__ ((malloc))
__attribute__ ((returns_nonnull))
__attribute__ ((assume_aligned(64)))
__attribute__ ((alloc_size(1)));
static inline
Ipp32f * gms_ippsMalloc_32f(int32_t len) {

    return (ippsMalloc_32f(len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__ ((malloc))
__attribute__ ((returns_nonnull))
__attribute__ ((assume_aligned(64)))
__attribute__ ((alloc_size(1)));
static inline
Ipp64f * gms_ippsMalloc_64f(int32_t len) {

    return (ippsMalloc_64f(len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__ ((malloc))
__attribute__ ((returns_nonnull))
__attribute__ ((assume_aligned(64)))
__attribute__ ((alloc_size(1)));
static inline
Ipp32fc * gms_ippsMalloc_32fc(int32_t len) {

     return (ippsMalloc_32fc(len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
__attribute__ ((malloc))
__attribute__ ((returns_nonnull))
__attribute__ ((assume_aligned(64)))
__attribute__ ((alloc_size(1)));
static inline
Ipp64fc * gms_ippsMalloc_64fc(int32_t len) {

     return (ippsMalloc_64fc(len));
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
static inline
void gms_ippsFree(void * ptr) {

     ippsFree(ptr);
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
static inline
bool gms_ippGetCacheParams() {
   /* Displays CPU cache information
      Based on official example.
    */
char* cacheType[] = {
    "Data Cache",
    "Instruction Cache",
    "Unified Cache"
};
   IPPcache * __restrict pCacheInfo;
   int32_t i;
   IPPStatus stat;
   stat = ippGetCacheParams(&pCacheInfo);
   if(stat != ippStsNoErr) {
      printf("Intel(R) Integrated Primitives (Intel(R) IPP) function returned error %s\n",
              ippGetStatusString( sts ));
      return (false);
   }
   i = 0;
   do{
        printf("cache type  = %s\n", cacheType[pCacheInfo[i].type-1] );
        printf("cache level = %d\n", pCacheInfo[i].level );
        printf("cache size  = %d\n", pCacheInfo[i].size );
        printf("+--------------------------------------+\n" );
    } while( pCacheInfo[++i].type > 0 );

    stat = ippGetL2CacheSize( &i );
    printf("\nCache L2 size = %d\n", i );
    return (true);
}


__ATTR_ALWAYS_INLINE__
__ATTR_HOT__
static inline
Ipp64u gms_ippGetCpuClocks() {

    return (ippGetCpuClocks());
}


__ATTR_ALWAYS_INLINE__
__ATTR_COLD__
static inline
IppStatus gms_ippGetCpuFreqMhz(int32_t * Mhz) {

    return (ippGetCpuFreqMhz(&Mhz));
}


#endif /*__GMS_IPPS_WRAPPERS_H__*/
