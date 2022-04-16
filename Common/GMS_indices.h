
#ifndef __GMS_INDICES_H__
#define __GMS_INDICES_H__

/*
	Flat multidiensional array indexing macros
*/

// To be used with CUDA
#ifndef INDEX2
#define INDEX2(isize,i,j) i + isize*j
#endif
#ifndef INDEX3
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#endif
#ifndef INDEX4
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#endif
#ifndef INDEX5
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))
#endif

/*
	@Warning:
				Macro parameters named idim,jdim,kdim,ldim
				must be present in the calling scope.
*/
//#define Idx2D(i,j)       ((i) * jdim + (j))

//#define Idx3D(i,j,k)     ((i) * jdim + ((j) * kdim + (k)))

//#define Idx4D(i,j,k,l)   ((i) + idim * ((j) + jdim * ((k) + kdim * (l))))

//#define Idx5D(i,j,k,l,m) ((i) + idim * ((j) + jdim * (k) + kdim * (l) + ldim * (m)))

/*
	@Warning:
				Macro parameters named nx,ny,nz must be
				present in the calling scope.
*/

//#define I2D(i,j)    ((i) * ny + (j))

//#define I3D(i,j,k)  ((i) * ny + ((j) * nz + (k)))

//#define I4D(i,j,k,l)  ((i) + nx * ((j) + ny * ((k) + nz * ((l)))))

/*
	In this case dimension parameters 
	must be passed to parametrized macro.
				
*/

#define Ix2D(i,ny,j) ((i) * (ny) + (j))

#define Ix3D(i,ny,j,nz,k) ((i) * (ny) + ((j) * (nz) + (k)))

// Different indexing scheme.


/*
	@Warning:
				Macro parameters named idim,jdim,kdim,ldim
				must be present in the calling scope.
*/
//#define IS2D(i,idim,j)   i * jdim + j

//#define IS3D(i,idim,j,jdim,k) i + idim * (j + jdim * k)

//#define IS4D(i,idim,j,jdim,k,kdim,l) i + idim * (j + jdim * (k + kdim * l))

//#define IS5D(i,idim,j,jdim,k,kdim,l,ldim,m) i + idim * (j + jdim * (k + kdim + (l + ldim * m)))

// WRF indexing scheme
// m_jme and m_kme must be present in the calling scope.
//#define Dim3(j,k,i) ((j) + (m_jme) * ((k) + (m_kme) * (i)))

/*
   int[] sixDArrayToOneDArray(final int[][][][][][] sixDArray) {

    int dim1Length = sixDArray.length;
    int dim2Length = sixDArray[0].length;
    int dim3Length = sixDArray[0][0].length;
    int dim4Length = sixDArray[0][0][0].length;
    int dim5Length = sixDArray[0][0][0][0].length;
    int dim6Length = sixDArray[0][0][0][0][0].length;

    int[] result = new int[dim1Length * dim2Length * dim3Length * dim4Length * dim5Length * dim6Length];

    for (int i1 = 0; i1 < dim1Length; i1++)
        for (int i2 = 0; i2 < dim2Length; i2++)
            for (int i3 = 0; i3 < dim3Length; i3++)
                for (int i4 = 0; i4 < dim4Length; i4++)
                    for (int i5 = 0; i5 < dim5Length; i5++)
                        for (int i6 = 0; i6 < dim6Length; i6++) {
                            int oneDIndex = i1 * dim2Length * dim3Length * dim4Length * dim5Length * dim6Length
                                          + i2 * dim3Length * dim4Length * dim5Length * dim6Length                                        
                                          + i3 * dim4Length * dim5Length * dim6Length
                                          + i4 * dim5Length * dim6Length
                                          + i5 * dim6Length
                                          + i6;
                            result[oneDIndex] = sixDArray[i1][i2][i3][i4][i5][i6];
                        }

    return result;
}
*/

#endif /*__GMS_INDICES_H__*/
