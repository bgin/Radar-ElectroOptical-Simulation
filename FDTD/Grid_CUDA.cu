/******************************************************************************
  * Copyright (c) 2023, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

#include "Grid_CUDA.hpp"
#include <cuda.h>
#include <cuda_runtime_api.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__device__ bool device_contain_point(double *vertx, double *verty, int vertlen, double ptx, double pty)
{
    bool inside = false;
    for(int i=0, j=vertlen-1; i< vertlen; j=i++){
        if( (verty[i]>pty)!=(verty[j]>pty) ){
            if(ptx < (vertx[j]-vertx[i])/(verty[j]-verty[i])*(pty-verty[i])+vertx[i]){
                inside=!inside;
            }
        }
    }
    return inside;
}


__device__ bool device_grid_intersect_poly(double *vertx, double *verty, double bbox0, double bbox1, double bbox2,
                            double bbox3, int vertlen, double xmin, double xmax, double ymin, double ymax)
{
    if(xmax<bbox0 || xmin>bbox1 || ymax<bbox2 || ymin>bbox3) { return false; }
    double xinter;
    for(int i=0, j=vertlen-1; i< vertlen; j=i++){
        if( (verty[i]>ymin)!=(verty[j]>ymin) )
        {
            xinter = (vertx[j]-vertx[i])/(verty[j]-verty[i])*(ymin-verty[i])+vertx[i]; if(xinter>=xmin && xinter<=xmax){ return true;}
        }
        if( (verty[i]>ymax)!=(verty[j]>ymax) )
        {
            xinter = (vertx[j]-vertx[i])/(verty[j]-verty[i])*(ymax-verty[i])+vertx[i]; if(xinter>=xmin && xinter<=xmax){ return true;}
        }
    }
    return false;
}

__device__ void device_overlap_raster_multicell(double *vertx, double *verty, int vertlen,
                            double xmin, double xminp1, double xmax, double ymin,
                             double yminp1, double ymax, bool *subcellflag, int Nsubcell, int Ncell, double *overlap,
                             int *inpoly)
{
    int Ninter=0; // number of intercepts
    int ub, lb;
    bool *vertin = new bool[vertlen];
    if(vertin==NULL){ 
        printf("bad_alloc vertin in device_overlap_raster_multicell\n"); 
        return; 
    }
    double *xinter = new double[vertlen];
    if(xinter==NULL){ 
        printf("bad_alloc xinter in device_overlap_raster_multicell\n"); 
        return; 
    }
    bool swapped;
    double subcellsize = xminp1-xmin, xtmp, pty;
    double subcellsizey = yminp1 - ymin;
    // detect if polygon vertex pairs span over the grid
    for(int i=0, j=vertlen-1; i< vertlen; j=i++){
        if( (verty[i]>ymax && verty[j]>ymax) || (verty[i]<ymin && verty[j]<ymin) ){ vertin[i]=false; }
        else{ vertin[i]=true; }
    }
    for(int i=0; i<Ncell; i++){ inpoly[i] = 0; }
    // horizontal ray tracing by y-subcell
    for(int cy=0; cy < Nsubcell; cy++){
        // rasterize grid subcell per each y to find out all x-intercepts
        pty = ymin + cy * subcellsizey;
        Ninter = 0;
        for(int i=0, j=vertlen-1; i< vertlen; j=i++){
            if(vertin[i]==true){
                if((verty[i]>pty)!=(verty[j]>pty)){
                    xinter[Ninter] = (vertx[j]-vertx[i])/(verty[j]-verty[i])*(pty-verty[i])+vertx[i];
                    Ninter++;
                }
            }
        }
        if(Ninter%2!=0){ printf("Odd number of interceptions.\n"); }
        if(Ninter>1){
            // bubble sort x-intercepts from left to right
            for(int k=1; k< Ninter; k++){
                swapped = false;
                for(int j=1; j< Ninter; j++){
                    if(xinter[j-1]>xinter[j]){xtmp=xinter[j-1]; xinter[j-1]=xinter[j]; xinter[j]=xtmp; swapped=true;}
                }
                if(swapped==false){goto ENDSORT;}
            }
            ENDSORT: ;
            // include grid subcell points in polygons
            for(int k=0; k< Ninter-1; k+=2){
                if(xinter[k] < xmin) { lb = 0;}
                else if(xinter[k] <= xmax) { lb = (int)ceil((xinter[k]-xmin)/subcellsize); }
                else if(xinter[k] > xmax) { goto ENDINCLUDE;}
                if(xinter[k+1] < xmin) {goto SKIPINTERVAL;}
                else if(xinter[k+1] < xmax) { ub = (int)ceil((xinter[k+1]-xmin)/subcellsize); }
                else if(xinter[k+1] >= xmax){ ub = Nsubcell*Ncell; }
                for(int cx=lb; cx<ub; cx++){
                    for(int ci=0; ci<Ncell; ci++){
                        if(cx < (ci+1)*Nsubcell && cx >= ci*Nsubcell){
                            if(subcellflag[cy*Ncell*Nsubcell+cx] == true){
                                inpoly[ci]++;
                                subcellflag[cy*Ncell*Nsubcell+cx] = false;
                            }
                        }
                    }
                }
                SKIPINTERVAL: ;
            }
            ENDINCLUDE: ;
        }
    }
    for(int i=0; i<Ncell; i++) { overlap[i] = 1.0*inpoly[i]/Nsubcell/Nsubcell; }
    delete[] xinter; delete[] vertin;
}


__global__ void kernel_get_values_Ncell(kernelpar *kpar)
{
    int Ncell;
    int i, j, k;
    k = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    i = blockIdx.z;
    int Nx0 = (int)ceil(kpar->Nx*1.0/kpar->NcellBlk);

//    if(k - kpar->k1 >= Nx0  || j - kpar->j1 >= kpar->Ny || i - kpar->i1 >= kpar->Nz) { return; } // non-laminated launch
    if(k >= Nx0  || j >= kpar->Ny || i >= kpar->Nz ) { return; }  // laminated launch
    if(k == Nx0 -1 && kpar->Nx%kpar->NcellBlk!=0) {
        if(kpar->Nx==1){ Ncell = 1; }
        else{ Ncell = kpar->Nx%kpar->NcellBlk; }
    }
    else{ Ncell = kpar->NcellBlk; }

    size_t ind_global = (kpar->i2 - kpar->i1) * kpar->Nx * kpar->Ny + j * kpar->Nx + k * kpar->NcellBlk;  // laminated launch

    if(ind_global >= kpar->size){ return;}

    double xmaxN;
    int Nsubcell = kpar->Nsubcell;
    double *cellarea = kpar->cellarea + ind_global;
    double *overlap = kpar->overlap + ind_global;
    bool *subcellflag = kpar->subcellflag + (size_t)j*kpar->Nx*Nsubcell*Nsubcell+k*kpar->NcellBlk*Nsubcell*Nsubcell;
    int *inpoly = kpar->inpoly + j * kpar->Nx + k * kpar->NcellBlk;
    double bbox0, bbox1, bbox2, bbox3;
    double cellmaxarea, cellarearatiomin, cellarearatiomax;
    double subcellx0, subcellx1, subcellxend, subcelly0, subcelly1, subcellyend;
    int nc;
    bool subcellflagset = false;
    int primidx=0;
    double zmin, zmax, ymin, ymax, xmin, xmax;
    double zratio = 1.0;
    bool containp1, containp2, containp3, containp4, gridcut;
    zmin = (kpar->i2+kpar->sz-0.5)*kpar->dz; zmax = (kpar->i2+kpar->sz+0.5)*kpar->dz;
    xmin = (k * kpar->NcellBlk + kpar->k1 + kpar->sx-0.5)*kpar->dx;
    xmax = (k * kpar->NcellBlk + kpar->k1 + kpar->sx+0.5)*kpar->dx;
    xmaxN = (k * kpar->NcellBlk + kpar->k1 + kpar->sx+0.5+(Ncell-1))*kpar->dx;
    ymin = ((j+kpar->j1)+kpar->sy-0.5)*kpar->dy; ymax = ((j+kpar->j1)+kpar->sy+0.5)*kpar->dy;

    cellmaxarea = (xmax-xmin)*(ymax-ymin);
    subcellx0 = xmin+(xmax-xmin)/Nsubcell*0.5;
    subcellx1 = xmin+(xmax-xmin)/Nsubcell*1.5;
    subcellxend = xmin+(xmax-xmin)/Nsubcell*(Ncell*Nsubcell-1+0.5);
    subcelly0 = ymin+(ymax-ymin)/Nsubcell*0.5;
    subcelly1 = ymin+(ymax-ymin)/Nsubcell*1.5;
    subcellyend = ymin+(ymax-ymin)/Nsubcell*(Nsubcell-1+0.5);
    if(zmax<=kpar->ZList[0]){
        for(int ci=0; ci<Ncell; ci++){ kpar->grid[ind_global+ci]=kpar->background; }
        printf("Simulation domain below Zmin.\n");
    }
    else{
        if(zmax>kpar->ZList[0] && zmin < kpar->ZList[0]){
            printf("zmax-zmin out-of-order\n");
            for(int ci=0; ci<Ncell; ci++) { kpar->grid[ind_global+ci] = (kpar->ZList[0]-zmin)/kpar->dz*1.0; }
            zmin = kpar->ZList[0];
        }
        for(int l=0; l<kpar->LayerLen; l++){
            if(zmin >= kpar->ZList[l] && zmax <= kpar->ZList[l+1]){
                zratio = (zmax-zmin)/kpar->dz;
                // criteria of using cache in previous grids:
                // 1. current grid does not contain interface
                // 2. not the first grid
                // 3. previous grid's starting Z is within the current material
                // 4. the previous grid's computation is finished.  Note this is a proxy block sync flag.  Index of
                //    the previous grid <0.5 implies that the computation in the previous grid is not finished.
                if(fabs(zratio-1.0)<1e-6 && (kpar->i2 > kpar->i1) && zmin-kpar->dz >= kpar->ZList[l]
                        && kpar->grid[ind_global-kpar->Nx*kpar->Ny].real()>0.5) {
                    for(int ci=0; ci<Ncell; ci++){
                        if(kpar->grid[ind_global+ci-kpar->Nx*kpar->Ny].real()>0.5){
                            kpar->grid[ind_global+ci] = kpar->grid[ind_global+ci-kpar->Nx*kpar->Ny];
                        }
                        else{ goto NOCACHE; }
                    }
                    goto ENDFORLOOP;
                } //use cache
                NOCACHE: ;

                for(int ci=0; ci<Ncell; ci++){ cellarea[ci] = cellmaxarea; }
                subcellflagset = false;
                primidx = kpar->PrimIdxPerLayer[l];
                for(int m=0; m<kpar->PrimNumberPerLayer[l]; m++){
                    bbox0 = kpar->PrimVertXmin[primidx]; if(bbox0>xmaxN){goto NOCOMP;}
                    bbox1 = kpar->PrimVertXmax[primidx]; if(bbox1<xmin){goto NOCOMP;}
                    bbox2 = kpar->PrimVertYmin[primidx]; if(bbox2>ymax){goto NOCOMP;}
                    bbox3 = kpar->PrimVertYmax[primidx]; if(bbox3<ymin){goto NOCOMP;}
                    // prim is rectangle AND grid falls within rectangle
                    if (subcellflagset==false && bbox0<=xmin && bbox1>=xmaxN && bbox2<=ymin && bbox3>=ymax &&
                            kpar->PrimVertNumber[primidx]==5){
                        for(int ci=0; ci<Ncell; ci++){
                            kpar->grid[ind_global+ci] += zratio * kpar->PrimMatValue[primidx];
                        }
                        goto ENDFORLOOP_RETURN;
                    }
                    if(xmin<bbox0 || xmin>bbox1 || ymin<bbox2 || ymin>bbox3){ containp1 = false;}
                    else{
                    containp1 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmin, ymin); }
                    if(xmaxN<bbox0 || xmaxN>bbox1 || ymin<bbox2 || ymin>bbox3){ containp2 = false;}
                    else{
                    containp2 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmaxN, ymin); }
                    if(xmaxN<bbox0 || xmaxN>bbox1 || ymax<bbox2 || ymax>bbox3){ containp3 = false;}
                    else{
                    containp3 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmaxN, ymax); }
                    if(xmin<bbox0 || xmin>bbox1 || ymax<bbox2 || ymax>bbox3){ containp4 = false;}
                    else{
                    containp4 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmin, ymax); }
                    gridcut = device_grid_intersect_poly(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        bbox0, bbox1, bbox2, bbox3,
                                        kpar->PrimVertNumber[primidx], xmin, xmaxN, ymin, ymax);
                    cellarearatiomin = 1.0;
                    for(int ci=0; ci<Ncell; ci++){
                        if(cellarea[ci]/cellmaxarea < cellarearatiomin){
                            cellarearatiomin = cellarea[ci]/cellmaxarea;
                        }
                    }
                    if(containp1 && containp2 && containp3 && containp4 && cellarearatiomin==1.0 && gridcut==false){
                        for(int ci=0; ci<Ncell; ci++) {
                            kpar->grid[ind_global+ci] += zratio * kpar->PrimMatValue[primidx];
                        }
                        goto ENDFORLOOP_RETURN;
                    }
                    else if( containp1 || containp2 || containp3 || containp4 || gridcut ){
                        if(subcellflagset==false){
                            memset(subcellflag, true, Ncell*Nsubcell*Nsubcell*sizeof(bool));
                            subcellflagset = true;
                        }
                        device_overlap_raster_multicell(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                            kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                            kpar->PrimVertNumber[primidx],
                                            subcellx0, subcellx1, subcellxend, subcelly0, subcelly1, subcellyend,
                                            subcellflag, Nsubcell, Ncell, overlap, inpoly);
                        for(int ci=0; ci<Ncell; ci++){
                            nc=0;
                            for(int cx=0; cx< Nsubcell; cx++){
                                for(int cy=0; cy<Nsubcell; cy++){
                                    if(subcellflag[cy*Ncell*Nsubcell+cx+ci*Nsubcell]==true){nc++;}
                                }
                            }
                            cellarea[ci] = cellmaxarea*nc/Nsubcell/Nsubcell;
                            kpar->grid[ind_global+ci] += overlap[ci] * zratio * kpar->PrimMatValue[primidx];
                        }
                    }
                    NOCOMP: ;
                    primidx++;
                    cellarearatiomax = 0.0;
                    for(int ci=0;ci<Ncell;ci++){
                        if(cellarea[ci]/cellmaxarea > cellarearatiomax){
                            cellarearatiomax=cellarea[ci]/cellmaxarea;
                        }
                    }
                    if(cellarearatiomax == 0) { goto ENDFORLOOP_RETURN;}
                } // for-m
                ENDFORLOOP_RETURN: ;
                goto ENDFORLOOP;
            } // if zmin
            else if(zmin >= kpar->ZList[l] && zmin < kpar->ZList[l+1] && zmax > kpar->ZList[l+1]){
                zratio = (kpar->ZList[l+1]-zmin) / kpar->dz;
                for(int ci=0; ci<Ncell; ci++){ cellarea[ci] = cellmaxarea; }
                subcellflagset = false;
                primidx = kpar->PrimIdxPerLayer[l];
                for(int m=0; m<kpar->PrimNumberPerLayer[l]; m++){
                    bbox0 = kpar->PrimVertXmin[primidx]; if(bbox0>xmaxN){goto NOCOMP2;}
                    bbox1 = kpar->PrimVertXmax[primidx]; if(bbox1<xmin){goto NOCOMP2;}
                    bbox2 = kpar->PrimVertYmin[primidx]; if(bbox2>ymax){goto NOCOMP2;}
                    bbox3 = kpar->PrimVertYmax[primidx]; if(bbox3<ymin){goto NOCOMP2;}
                    // prim is rectangle AND grid falls within rectangle
                    if (subcellflagset==false && bbox0<=xmin && bbox1>=xmaxN && bbox2<=ymin && bbox3>=ymax &&
                            kpar->PrimVertNumber[primidx]==5){
                        for(int ci=0; ci<Ncell; ci++){
                            kpar->grid[ind_global+ci] += zratio * kpar->PrimMatValue[primidx];
                        }
                        goto ENDFORLOOP_RETURN2;
                    }
                    if(xmin<bbox0 || xmin>bbox1 || ymin<bbox2 || ymin>bbox3){ containp1=false;}
                    else{
                    containp1 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmin, ymin);}
                    if(xmaxN<bbox0 || xmaxN>bbox1 || ymin<bbox2 || ymin>bbox3){ containp2=false;}
                    else{
                    containp2 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmaxN, ymin);}
                    if(xmaxN<bbox0 || xmaxN>bbox1 || ymax<bbox2 || ymax>bbox3){ containp3=false;}
                    else{
                    containp3 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmaxN, ymax);}
                    if(xmin<bbox0 || xmin>bbox1 || ymax<bbox2 || ymax>bbox3){ containp4=false;}
                    else{
                    containp4 = device_contain_point(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertNumber[primidx], xmin, ymax);}
                    gridcut = device_grid_intersect_poly(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                        kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                        bbox0, bbox1, bbox2, bbox3,
                                        kpar->PrimVertNumber[primidx], xmin, xmaxN, ymin, ymax);
                    cellarearatiomin = 1.0;
                    for(int ci=0; ci<Ncell; ci++){
                        if(cellarea[ci]/cellmaxarea < cellarearatiomin){
                            cellarearatiomin = cellarea[ci]/cellmaxarea;
                        }
                    }
                    if(containp1 && containp2 && containp3 && containp4 && cellarearatiomin==1.0 && gridcut==false){
                        for(int ci=0; ci<Ncell; ci++){
                            kpar->grid[ind_global+ci] += zratio * kpar->PrimMatValue[primidx];
                        }
                        goto ENDFORLOOP_RETURN2;
                    }
                    else if(containp1 || containp2 || containp3 || containp4 || gridcut){
                        if(subcellflagset==false){
                            memset(subcellflag, true, Ncell*Nsubcell*Nsubcell*sizeof(bool));
                            subcellflagset=true;
                        }
                        device_overlap_raster_multicell(kpar->PrimVertXAll+kpar->PrimVertCoordLoc[primidx],
                                            kpar->PrimVertYAll+kpar->PrimVertCoordLoc[primidx],
                                            kpar->PrimVertNumber[primidx],
                                            subcellx0, subcellx1, subcellxend, subcelly0, subcelly1, subcellyend,
                                            subcellflag, Nsubcell, Ncell, overlap, inpoly);
                        for(int ci=0; ci<Ncell; ci++){
                            nc=0;
                            for(int cx=0; cx< Nsubcell; cx++){
                                for(int cy=0; cy<Nsubcell; cy++){
                                    if(subcellflag[cy*Ncell*Nsubcell+cx+ci*Nsubcell]==true){nc++;}
                                }
                            }
                            cellarea[ci] = cellmaxarea*nc/Nsubcell/Nsubcell;
                            kpar->grid[ind_global+ci] += zratio * overlap[ci] * kpar->PrimMatValue[primidx];
                        }
                    }
                    NOCOMP2: ;
                    primidx++;
                    cellarearatiomax = 0.0;
                    for(int ci=0;ci<Ncell;ci++){
                        if(cellarea[ci]/cellmaxarea>0){
                            cellarearatiomax=cellarea[ci]/cellmaxarea;
                        }
                    }
                    if(cellarearatiomax==0) { goto ENDFORLOOP_RETURN2;}
                }  // for m
                ENDFORLOOP_RETURN2: ;
                zmin = kpar->ZList[l+1];
            }  // else if
        }  // for l
    }  // end else
    ENDFORLOOP: ;
}

void invoke_kernel_get_values_Ncell(kernelpar *kpar, dim3 grid_dim, dim3 block_dim)
{
    kernel_get_values_Ncell <<< grid_dim, block_dim >>> (kpar);
}
