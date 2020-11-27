
#ifndef __ARMA_ESTIMATION_HPP__
#define __ARMA_ESTIMATION_HPP__


#include <cstdint>
#include <cmath>

/*
       Converted from F90 implementation of the following ASA algorithms
      ALGORITHM AS 154  APPL. STATIST. (1980) VOL.29, P.311
      ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2
*/

/*
   !  ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2

!  Finite sample prediction from ARIMA processes.
*/

__attribute__((aligned(32)))
__attribute__((always_inline))
__attribute__((hot))
static inline
void forkal(const int32_t ip,
            const int32_t iq,
            const int32_t ir,
            const int32_t np,
            const int32_t ird,
            const int32_t irz,
            const int32_t id,
            const int32_t il,
            const int32_t n,
            const int32_t nrbar,
            float * __restrict __attribute__((aligned(64))) phi,
            const float * __restrict __attribute__((aligned(64))) theta,
            const float * __restrict __attribute__((aligned(64))) delta,
            float * __restrict __attribute__((aligned(64))) w,
            float * __restrict __attribute__((aligned(64))) y,
            float * __restrict __attribute__((aligned(64))) amse,
            float * __restrict __attribute__((aligned(64))) a,
            float * __restrict __attribute__((aligned(64))) p,
            float * __restrict __attribute__((aligned(64))) v,
            float * __restrict __attribute__((aligned(64))) resid,
            float * __restrict __attribute__((aligned(64))) e,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) store,
            int32_t &ifault ) {
      //
    // Invoking this routine will calculate the finite sample predictions
     //     and their conditional mean square errors for any ARIMA process.
     // Will surely impact the branch predictor resources.
     int32_t k;
     ifault = 0;
     if(ip<0) ifault = 1;
     if(iq<0) ifault = ifault + 2;
     if(ip*ip+iq*iq == 0) ifault = 4;
     k = iq + 1;
     if(k<ip) k = ip;
     if(ir!=k) ifault = 5;
     if(np!=ir*(ir+1)/2) ifault = 6;
     if(nrbar!=np*(np-1)/2) ifault = 7;
     if(id<0) ifault = 8;
     if(ird!=ir+id) ifault = 9;
     if(irz!=ird*(ird+1)/2) ifault = 10;
     if(il<1) ifault = 11;
     if(ifault!=0) return;

     int32_t i, i45, ibc, id1, id2r, id2r1, id2r2, idd1, idd2, iddr, idrr1, 
            idk, iid, ind, ind1, ind2, iq1, ir1, ir2, iri, iri1, irj, iupd, 
            j, j1, jj, jkl, jkl1, jklj, jrj, jrk, k, k1, kk, kk1, kkk, l,
            lk, lk1, ll, lli, nit, nj, nt;
     float a1, aa, ams, del, dt, phii, phij, phijdt, sigma, ssq, sumlog;
     //  Calculate initial conditions for Kalman filter
     __assume_aligned(a,64);
     __assume_aligned(v,64);
     a[0] = 0.0f;
     v[0] = 0.0f;
     if(np==1) goto L130;
#pragma vector aligned
#pragma vector always
     for(int32_t idx = 2; idx != np; ++idx) {
         v[i] = 0.0f;
     }
     if(iq==0) goto L130;
     iq1 = iq+1;
     __assume_aligned(theta,64);
#pragma ivdep
#pragma vector aligned
#pragma vector always
     for(i = 1; i != iq1; ++i) {
         v[i] = theta[i-1];
     }
     for(j = 0; j != iq; ++j) {
        ll = j * (ir+ir+1-j) / 2;
        float tj = theta[j];
        for(i = j; i != iq; ++i) {
            lli = ll + i;
            v[lli] = theta[i] * tj;
        }
    }
    /*
         Find initial likelihood conditions.
!     IFAULT not tested on exit from STARMA as all possible errors
!     have been checked above.
    */
L130:
     if(ir==1) p[1] = 1.0f/(1.0f-phi[0]*phi[0]);
     if(ir!=1) {
         starma(ip,iq,ir,np,phi,theta,a,p,v,
                thetab,xnext,xrow,rbar,nrbar,ifault);
     }  
     //    Calculate data transformations
     nt = n - id;
     if(id==0) goto L170;
     __assume_aligned(store,64);
     __assume_aligned(w,64);
#pragma vector aligned
#pragma vector always
     for(j = 1; j != id; ++j) {
         nj = n - j;
         store[j] = w[nj];
     }
     __assume_aligned(delta,64);
     for(i = 0; i != nt; ++i) {
         aa = 0.0f;
         for(k = 0; k != id; ++k) {
             idk = id + i - k;
             aa = aa - delta[k] * w[idk];
         }
         iid = i + id;
         w[i] = w[iid] + aa;
     }
     //  Evaluate likelihood to obtain final KF conditions
L170:
     sumlog = 0.0f;
     ssq    = 0.0f;
     iupd   = 1;
     del    = -1.0f;
     nit    = 0;
     karma(ip,iq,ir,phi,theta,a,p,v,nt,w,resid,
           sumlog,ssq,iupd,del,e,nit);
     //  Calculate M.L.E. of sigma squared
     sigma = 0.0f;
     __assume_aligned(resid,64);
#pragma vector aligned
#pragma simd reduction(+:sigma)
     for(j = 0; j != nt; ++j) {
         sigma = sigma+resid[j]*resid[j];
     }
     // Reset the initial A and P when differencing occurs
     if(id==0) goto L250;
     __assume_aligned(xrow,64);
     __assume_aligned(p,64);
#pragma vector aligned
#pragma vector always
     for(i = 0; i != np; ++i) {
         xrow[i] = p[i];
     }   
#pragma vector aligned
#pragma vector always
     for(i = 0; i != irz; ++i) {
         p[i] = 0.0f;
     } 
     ind = 0;
     
}

/*
    !  INVOKING THIS SUBROUTINE SETS THE VALUES OF V AND PHI, AND
!  OBTAINS THE INITIAL VALUES OF A AND P.
!  THIS ROUTINE IS NOT SUITABLE FOR USE WITH AN AR(1) PROCESS.
!  IN THIS CASE THE FOLLOWING INSTRUCTIONS SHOULD BE USED FOR INITIALISATION.
!  V(1) = 1.0
!  A(1) = 0.0
!  P(1) = 1.0 / (1.0 - PHI(1) * PHI(1))
*/

__attribute__((aligned(32)))
__attribute__((always_inline))
__attribute__((hot))
static inline
void starma(const int32_t ip,
            const int32_t iq,
            const int32_t ir,
            const int32_t np,
            float * __restrict __attribute__((aligned(64))) phi,
            const float * __restrict attribute__((aligned(64))) theta,
            float * __restrict __attribute__((aligned(64))) a,
            float * __restrict __attribute__((aligned(64))) p,
            float * __restrict __attribute__((aligned(64))) v,
            float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            float * __restrict __attribute__((aligned(64))) rbar,
            const int32_t nrbar,
            int32_t &ifault) {
    //  CHECK FOR FAILURE INDICATION.
    // The chain of 10 branches likely a performance hit.
    // THere should not be more than 15 branches in chain of
    // instructions.
    // Forward branches predicted not-taken.
    int32_t k;
    ifault = 0;
    if(ip<0) { ifault = 1;} // likely not-taken
    if(iq<0) { ifault += 2;} // likely not-taken
    if(ip==0 && iq==0) { ifault = 4;} // likely not-taken
    k = iq + 1;
    if(k<ip) { k = ip;}
    if(ir!=k) { ifault = 5;}
    if(np!=ir*(ir+1)/2) { ifault = 6;}
    if(nrbar!=np*(np-1)/2) { ifault = 7;}
    if(ir==1) { ifault = 8;}
    if(ifault!=0) { return;}
    //loat vj,phii,phij,ssqerr,recres,ynext;
    register float phii,phij,ynext;
    float ssqer,recres;
    int32_t ifail,irank,npr,npr1
    register int32_t ind,ind1,ind2,indj,indi;
    // NOW SET A(0), V AND PHI.
    __assume_aligned(phi,64);
    __assume_aligned(theta,64);
    __assume_aligned(a,64);
    __assume_aligned(p,64);
    __assume_aligned(v,64);
    __assume_aligned(thetab,64);
    __assume_aligned(xnex,64);
    __assume_aligned(xrow,64);
    __assume_aligned(rbar,64);
#pragma ivdep
#pragma vector aligned
#pragma vector always
     for(int32_t i = 1; i != ir; ++i) {
         a[i] = 0.0f;
         if(i>ip) phi[i] = 0.0f;
         v[i] = 0.0f;
         if(i<=iq+1) v[i] = theta[i-1];
     }   
     a[0] = 0.0f;
     if(ip==0) phi[0] = 0.0f;
     v[0] = 0.0f;
     ind = ir;
     for(int32_t j = 1; j != ir; ++j) {
         vj = v[j];

#pragma vector aligned
#pragma vector always
        for(int32_t i = j; i != ir; ++i) {
            ind += 1;
            v[iind] = v[i] * vj;
        }
     }

    if(ip!=0) {
     //  NOW FIND P(0).
     /*
    
         THE SET OF EQUATIONS S * VEC(P(0)) = VEC(V) IS SOLVED FOR VEC(P(0)).
!   S IS GENERATED ROW BY ROW IN THE ARRAY XNEXT.
!   THE ORDER OF ELEMENTS IN P IS CHANGED, SO AS TO
!   BRING MORE LEADING ZEROS INTO THE ROWS OF S,
!   HENCE ACHIEVING A REDUCTION OF COMPUTING TIME.
     */
          irank = 0;
          ssqer = 0.0f;
#pragma vector aligned
#pragma vector always
        for(int32_t i = 0; i != nrbar; ++i) { rbar[i] = 0.0f;}
#pragma vector aligned
#pragma vector always  
            for(int32_t i = 0; i != np; ++i) {
                p[i] = 0.0f;
                thetab[i] = 0.0f;
                xnext[i] = 0.0f;
            }
            ind  = 0;
            ind1 = 0;
            npr  = np - ir;
            npr1 = npr + 1;
            indj = npr1;
            ind2 = npr;
            for(int32_t j = 0; j != ir; ++j) {
                phij = phi[j];
                xnext[indj] = 0.0f;
                indj += 1;
                indi = npr1 + j;
                for(int32_t i = j; i != ir; ++i) {
                     ind += 1;
                     ynext = v[ind];
                     phii = phi[i];
                     if(j != ir) {
                        xnext[indj] = -phii;
                        if(i != ir) {
                           xnexr[indi] = xnext[indi] - phij;
                           ind1 += 1;
                           xnext[ind1] = -1.0f;
                        }
                 }
                 xnext[npr1] = -phii * phij;
                 ind2 = ind2 + 1;
                 if(ind2>np) ind2 = 1;
                 xnext[ind2] = xnext[ind2] + 1.0f;
                 inclu2(np,1.0f,xnext,xrow,ynext,p,rbar,thetab,
                   ssqer,recres,irank,ifault);
                 // NO NEED TO CHECK IFAIL AS WEIGHT = 1.0
                 xnext[ind2] = 0.0f;
                 if(i != ir) {
                    xnext[indi] = 0.0f;
                    indi += 1;
                    xnext[ind1] = 0.0f;
                }
            }
        }
          regres(np,nrbar,rbar,thetab,p);
        // NOW RE-ORDER P.
          ind = npr;
#pragma vector aligned
#pragma vector always  
          for(int32_t i = 0; i != ir; ++i) {
              ind += 1;
              xnex[i] = p[ind];
           }
           ind = np;
           ind1 = npr;
           for(int3_t i = 0; i != npr; ++i) {
               p[ind] = p[ind1];
               ind  -= 1;
               ind1 -= 1;
           }
#pragma vector aligned
#pragma vector always 
           for(int32_t i = 0; i != ir; ++i) {
               p[i] = xnex[i];
           } 
           return;
       }
       //  P(0) IS OBTAINED BY BACKSUBSTITUTION FOR A MOVING AVERAGE PROCESS.
       indn = np+1;
       ind  = np+1;
       for(int32_t i = 0; i != ir; ++i) {
           for(int32_t j = 0; j != i; ++j) {
               ind -= 1;
               p[ind] = v[ind];
               if(j != 1) {
                   indn -= 1;
                   p[ind] += p[indn];
               }
           }
       }
}

/*
   ! N.B. Argument NRBAR has been removed.

!   ALGORITHM AS 154.3  APPL. STATIST. (1980) VOL.29, P.311

!   FORTRAN VERSION OF REVISED VERSION OF ALGORITHM AS 75.1
!   APPL. STATIST. (1974) VOL.23, P.448
!   SEE REMARK AS R17 APPL. STATIST. (1976) VOL.25, P.323
*/

__attribute__((aligned(32)))
__attribute__((always_inline))
__attribute__((hot))
static inline
void inclu2(const int32_t np,
            const float weight,
            const float * __restrict __attribute__((aligned(64))) xnext,
            float * __restrict __attribute__((aligned(64))) xrow,
            const float ynext,
            float * __restrict __attribute__((aligned(64))) d,
            float * __restrict __attribute__((aligned(64))) rbar,
            float * __restrict __attribute__((aligned(64))) thetab,
            float &ssqerr,
            float &recres,
            int32_t &irank,
            int32_t &ifault) {
    /*
       !   INVOKING THIS SUBROUTINE UPDATES D, RBAR, THETAB, SSQERR
!   AND IRANK BY THE INCLUSION OF XNEXT AND YNEXT WITH A
!   SPECIFIED WEIGHT. THE VALUES OF XNEXT, YNEXT AND WEIGHT WILL
!   BE CONSERVED.  THE CORRESPONDING VALUE OF RECRES IS CALCULATED.

    */
    ifault = 1;
    if(weigth<=0.0f) { return;}
    register float xk,rbthis,xk,xi,cbar,sbar;
    float wt,y,di,dpi;
    register int32_t ithisr,i1;
    y = ynext;
    wt = weight;
    if(wt<=0.0f)
    __assume_aligned(xrow,64);
    __assume_aligned(xnext,64);
#pragma vector aligned
#pragma vector always
     for(int32_t i = 0; i != np; ++i) {
         xrow[i] = xnext[i];
     } 
    recres = 0.0f;
    ifault = 0;
    ithisr = 0;
    __assume_aligned(d,64);
    __assume_aligned(rbar,64);
    __assume_aligned(thetab,64);
    for(int32_t i = 0; i != np; ++i) {
        if(xrow[i]==0.0f) {
            ithisr = ithisr + np - 1;
        }
        else {
            xi = xrow[i];
            di = d[i];
            dpi = di + wt * xi * xi;
            d[i] = dpi;
            cbar = di / dpi;
            sbar = wt * xi / dpi;
            wt = cbar * wt;
            if(i != np) {
                i1 = i + 1;
#pragma vector aligned
#pragma vector always
                for(int32_t k = i1; k != np; ++i1) {
                    ithisr += 1;
                    xk = xrow[k];
                    rbthis = rbar[ithisr];
                    xrow[k] = xk - xi * rbthis;
                    rbar[ithisr] = cbar * rbthis + sbar * xk;
                }
            }
            xt = y;
            y = xk - xi * thetab[i];
            thetab[i] = cbar * thetab[i] + sbar * xk;
            if(di==0.0f) goto lab_40;
        }
    }
    ssqerr = ssqerr + wt * y * y;
    recres = y * std::sqrt(wt);
    return;
lab_40:
    irank += 1;
}

/*
   !   ALGORITHM AS 154.4  APPL. STATIST. (1980) VOL.29, P.311

!   REVISED VERSION OF ALGORITHM AS 75.4
!   APPL. STATIST. (1974) VOL.23, P.448
!   INVOKING THIS SUBROUTINE OBTAINS BETA BY BACKSUBSTITUTION
!   IN THE TRIANGULAR SYSTEM RBAR AND THETAB.
*/

__attribute__((aligned(32)))
__attribute__((always_inline))
__attribute__((hot))
static inline
void regres(const int32_t np,
            const int32_t nrbar,
            const float * __restrict __attribute__((aligned(64))) rbar,
            const float * __restrict __attribute__((aligned(64))) thetab,
            float * __restrict __attribute__((aligned(64))) beta) {
    register float bi;
    register int32_t im,jm,ithisr;
    __assume_aligned(rbar,64);
    __assume_aligned(thetab,64);
    __assume_aligned(beta,64);
    ithisr = nrbar;
    im = np;
    for(int32_t i = 0; i != np; ++i) {
        bi = thetab[im];
        if(im!=np) {
            i1 = i - 1;
            jm = np;
#pragma vector aligned
#pragma vector always
            for(int32_t j = 0; j != i1; ++j) {
                bi = bi - rbar[ithisr] * beta[jm];
                ithisr -= 1;
                jm -= 1;
            }
        }
        beta[im] = bi;
        im -= 1;
    }
}









#endif /*__ARMA_ESTIMATION_HPP__*/