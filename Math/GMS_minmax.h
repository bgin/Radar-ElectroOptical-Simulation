
#ifndef __GMS_MINMAX_H__
#define __GMS_MINMAX_H__


#if defined(__GNUC__) || defined(__clang__)
#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})
#else
#error "Unsupported statement expression!!"
#endif











#endif /*__GMS_MINMAX_H__*/
