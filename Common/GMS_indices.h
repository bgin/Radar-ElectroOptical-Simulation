
#ifndef __GMS_INDICES_H__
#define __GMS_INDICES_H__

/*
	Flat multidiensional array indexing macros
*/

/*
	@Warning:
				Macro parameters named idim,jdim,kdim,ldim
				must be present in the calling scope.
*/
#define Idx2D(i,j)       ((i) * jdim + (j))

#define Idx3D(i,j,k)     ((i) * jdim + ((j) * kdim + (k)))

#define Idx4D(i,j,k,l)   ((i) + idim * ((j) + jdim * ((k) + kdim * (l))))

#define Idx5D(i,j,k,l,m) ((i) + idim * ((j) + jdim * (k) + kdim * (l) + ldim * (m)))

/*
	@Warning:
				Macro parameters named nx,ny,nz must be
				present in the calling scope.
*/

#define I2D(i,j)    ((i) * ny + (j))

#define I3D(i,j,k)  ((i) * ny + ((j) * nz + (k)))

#define I4D(i,j,k,l)  ((i) + nx * ((j) + ny * ((k) + nz * ((l)))))

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
#define IS2D(i,idim,j)   i * jdim + j

#define IS3D(i,idim,j,jdim,k) i + idim * (j + jdim * k)

#define IS4D(i,idim,j,jdim,k,kdim,l) i + idim * (j + jdim * (k + kdim * l))

#define IS5D(i,idim,j,jdim,k,kdim,l,ldim,m) i + idim * (j + jdim * (k + kdim + (l + ldim * m)))

// WRF indexing scheme
// m_jme and m_kme must be present in the calling scope.
#define Dim3(j,k,i) ((j) + (m_jme) * ((k) + (m_kme) * (i)))



#endif /*__GMS_INDICES_H__*/
