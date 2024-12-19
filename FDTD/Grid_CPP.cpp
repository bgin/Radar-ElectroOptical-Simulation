/******************************************************************************
 * Copyright (c) 2023, Andrew Michaels.  All rights reserved.
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

#include "Grid_CPP.hpp"
#include "Grid_CUDA.hpp"
#include <iostream>
#include <climits>
#include <ctime>
#include <exception>
#include <omp.h>

#include <stdio.h>

#include <fstream>
#include <iomanip>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
using namespace Grid;

/**************************************** Materials ****************************************/
//------------------------------ Grid Material ------------------------------------/
GridMaterial2D::GridMaterial2D(int M, int N, ArrayXXcd grid) : _M(M), _N(N), _grid(grid) {}

std::complex<double> GridMaterial2D::get_value(double x, double y)
{
    int xx = int(x), yy = int(y);
	if(xx > 0 && xx < _M && yy > 0 && yy < _N)
		return _grid(yy,xx);
	else
		return 1.0;
}

void GridMaterial2D::get_values(ArrayXcd& grid, int k1, int k2, int j1, int j2, double sx, double sy)
{
    int N = k2 - k1;
    for(int i = j1; i < j2; i++) {
        for(int j = k1; j < k2; j++) {
            grid((i-j1)*N + j-k1) = _grid(i,j);
        }
    }
}

void GridMaterial2D::set_grid(int M, int N, ArrayXXcd grid)
{
	_M = M;	_N = N;	_grid = grid;
}

int GridMaterial2D::get_M() { return _M; }
int GridMaterial2D::get_N() { return _N; }

//------------------------------ MaterialPrimitives ------------------------------------/
GridCell::GridCell()
{
}
		
void GridCell::set_vertices(double xmin, double xmax, double ymin, double ymax)
{
	double area;
	_verts.clear();

	Polygon_2D new_poly;
	boost::geometry::append(new_poly, Point_2D(xmin, ymin));
	boost::geometry::append(new_poly, Point_2D(xmin, ymax));
	boost::geometry::append(new_poly, Point_2D(xmax, ymax));
	boost::geometry::append(new_poly, Point_2D(xmax, ymin));
	boost::geometry::append(new_poly, Point_2D(xmin, ymin));
    boost::geometry::correct(new_poly);

    boost::geometry::assign(_original, new_poly);

    _verts.push_back(new_poly);
	area = boost::geometry::area(new_poly);
	_area = fabs(area);
	_max_area = _area;
}

double GridCell::intersect(const Polygon_2D poly)
{
	double area = 0.0,
		   intersected_area,
           geo_area;
	_diffs.clear();
	std::vector<Polygon_2D>::const_iterator i;

    // Do the difference
	for(i = _verts.begin(); i != _verts.end(); ++i) {
		boost::geometry::difference((*i), poly, _diffs);
	}
	_verts.clear();
	
	for(i = _diffs.begin(); i != _diffs.end(); i++) {
		_verts.push_back(*i);
        geo_area = boost::geometry::area(*i);
		area += fabs(geo_area);
	}

	intersected_area = _area - area;
	_area = area;
	return intersected_area/_max_area;
}

double GridCell::get_area()
{
	return _area;
}

double GridCell::get_max_area()
{
	return _max_area;
}

double GridCell::get_area_ratio()
{
	return _area/_max_area;
}		

//--------------------------------------------------------------------------------------/
//------------------------------ MaterialPrimitives ------------------------------------/
//--------------------------------------------------------------------------------------/
MaterialPrimitive::MaterialPrimitive()
{
	_layer = 1;
}

int MaterialPrimitive::get_layer() const { return _layer; }
void MaterialPrimitive::set_layer(int layer) { _layer = layer; }

bool MaterialPrimitive::operator<(const MaterialPrimitive& rhs)
{
	return _layer < rhs.get_layer();
}

Circle::Circle(double x0, double y0, double r) : _x0(x0), _y0(y0), _r(r) {}
Circle::~Circle() {}

bool Circle::contains_point(double x, double y)
{
	double dx = x - _x0,
		   dy = y - _y0;
	return dx*dx + dy*dy < _r*_r;
}

bool Circle::bbox_contains_point(double x, double y)
{
	double xmin = _x0-_r,
		   xmax = _x0+_r,
		   ymin = _y0-_r,
		   ymax = _y0+_r;

	return (x > xmin) && (x < xmax) && (y > ymin) && (y < ymax);
}

std::complex<double> Circle::get_material(double x, double y)
{
	return _mat;
}

// TODO: Implement this properly
double Circle::get_cell_overlap(GridCell& cell)
{
	return 0.0;
}	

void Circle::get_vertices(double *vx, double *vy)
{
    int p=1;
//     double *x, *y;
//     x = (double *)malloc(_verts.outer().size() * sizeof(double));
//     y = (double *)malloc(_verts.outer().size() * sizeof(double));
//     int idx = 0;
//     for(const auto& point : _verts.outer()){
//         x[idx] = point.x();
//         y[idx] = point.y();
//         idx++;
//     }
//     std::memcpy(dx, x, _verts.outer().size());
//     std::memcpy(dx, x, _verts.outer().size());
//     delete[] x; delete[] y;
}

int Circle::get_number_vertices()
{
    return 37;
}

void Circle::set_material(std::complex<double> mat)
{
	_mat = mat;
}

void Circle::set_position(double x0, double y0)
{
	_x0 = x0;
	_y0 = y0;
}

void Circle::set_radius(double r)
{
	_r = r;
}

double Circle::get_x0()
{
	return _x0;
}

double Circle::get_y0()
{
	return _y0;
}

double Circle::get_r()
{
	return _r;
}

Rectangle::Rectangle(double x0, double y0, double width, double height) : 
		_x0(x0), _y0(y0), _width(width), _height(height)
{
    // Points must be defined clockwise and the polygon must be closed
    boost::geometry::append(_verts, Point_2D(x0-width/2.0, y0-height/2.0));
    boost::geometry::append(_verts, Point_2D(x0-width/2.0, y0+height/2.0));
    boost::geometry::append(_verts, Point_2D(x0+width/2.0, y0+height/2.0));
    boost::geometry::append(_verts, Point_2D(x0+width/2.0, y0-height/2.0));
    boost::geometry::append(_verts, Point_2D(x0-width/2.0, y0-height/2.0));
}

Rectangle::~Rectangle()
{
}

bool Rectangle::contains_point(double x, double y)
{
	double hwidth = _width/2,
		   hheight = _height/2;
	return (x > _x0-hwidth) && (x < _x0+hwidth) && (y > _y0-hheight) && (y < _y0+hheight);
}

bool Rectangle::bbox_contains_point(double x, double y)
{
	return contains_point(x,y);
}

std::complex<double> Rectangle::get_material(double x, double y)
{
	return _mat;
}

double Rectangle::get_cell_overlap(GridCell& cell)
{
	return cell.intersect(_verts);
}

void Rectangle::get_vertices(double *vx, double *vy)
{
    int idx = 0;
    for(const auto& point : _verts.outer()){
        vx[idx] = point.x();
        vy[idx] = point.y();
        idx++;
    }
}

int Rectangle::get_number_vertices()
{
    return _verts.outer().size();
}

void Rectangle::set_material(std::complex<double> mat)
{
	_mat = mat;
}

void Rectangle::set_width(double w)
{
	_width = w;
    std::vector<Point_2D>& outer = _verts.outer();
    outer[0].x(_x0-_width/2.0);
    outer[1].x(_x0-_width/2.0);
    outer[2].x(_x0+_width/2.0);
    outer[3].x(_x0+_width/2.0);
    outer[4].x(_x0-_width/2.0);
}

void Rectangle::set_height(double h)
{
	_height = h;
    std::vector<Point_2D>& outer = _verts.outer();
    outer[0].y(_y0-_height/2.0);
    outer[1].y(_y0+_height/2.0);
    outer[2].y(_y0+_height/2.0);
    outer[3].y(_y0-_height/2.0);
    outer[4].y(_y0-_height/2.0);
}

void Rectangle::set_position(double x0, double y0)
{
	_x0 = x0;
	_y0 = y0;

    std::vector<Point_2D>& outer = _verts.outer();
    outer[0].x(_x0-_width/2.0);
    outer[1].x(_x0-_width/2.0);
    outer[2].x(_x0+_width/2.0);
    outer[3].x(_x0+_width/2.0);
    outer[4].x(_x0-_width/2.0);

    outer[0].y(_y0-_height/2.0);
    outer[1].y(_y0+_height/2.0);
    outer[2].y(_y0+_height/2.0);
    outer[3].y(_y0-_height/2.0);
    outer[4].y(_y0-_height/2.0);
}

//------------------------------ Polygon ------------------------------------/

Polygon::Polygon()
{
}

Polygon::Polygon(double* x, double* y, int n)
{
	set_points(x, y, n);
}

Polygon::~Polygon()
{
	_verts.clear();
}

void Polygon::add_point(double x, double y)
{
	boost::geometry::append(_verts, boost::geometry::make<Point_2D>(x,y));
    // update the bounding box
    boost::geometry::envelope(_verts, _bbox);
    // correct the geometry
    //boost::geometry::correct(_verts);
}

/**
 * NOTE: Currently a copy of the input points is made.  This will be slowish.
 */
void Polygon::add_points(double* x, double* y, int n)
{
    for(int i = 0; i < n; i++) {
        boost::geometry::append(_verts, boost::geometry::make<Point_2D>(x[i], y[i]));
    }
    // update the bounding box
    boost::geometry::envelope(_verts, _bbox);
    // correct the geometry
    boost::geometry::correct(_verts);
}

void Polygon::set_point(double x, double y, int index)
{
	Point_2D& p = _verts.outer()[index];
    p.x(x);
    p.y(y);
    // update the bounding box
    boost::geometry::envelope(_verts, _bbox);
    // assume the geometry is correct.
}

void Polygon::set_points(double* x, double* y, int n)
{
	_verts.clear();
	add_points(x,y,n);
    // update the bounding box
    boost::geometry::envelope(_verts, _bbox);
    // correct the geometry
    boost::geometry::correct(_verts);
}

bool Polygon::contains_point(double x, double y)
{
    Point_2D p(x, y);
    bool inside = false;

    if(bbox_contains_point(x,y)){
        if(boost::geometry::within(p, _verts)){
            inside = true;
        }
    }
	return inside;
}

bool Polygon::bbox_contains_point(double x, double y)
{
    Point_2D_D p(x,y);
    bool inside = boost::geometry::within(p, _bbox);
	return inside;
}

std::complex<double> Polygon::get_material(double x, double y)
{
	//if(contains_point(x,y))
	return _mat;
}

double Polygon::get_cell_overlap(GridCell& cell)
{
	return cell.intersect(_verts);
}

void Polygon::get_vertices(double *vx, double *vy)
{
    int idx = 0;
    for(const auto& point : _verts.outer()){
        vx[idx] = point.x();
        vy[idx] = point.y();
        idx++;
    }
}

int Polygon::get_number_vertices()
{
    return _verts.outer().size();
}

void Polygon::set_material(std::complex<double> mat)
{
	_mat = mat;
}

/////////////////////////////////////////////////////////////////////////////////////
// StructuredMaterial2D
/////////////////////////////////////////////////////////////////////////////////////
StructuredMaterial2D::StructuredMaterial2D(double w, double h, double dx, double dy) :
	_w(w), _h(h), _dx(dx), _dy(dy)
{}

StructuredMaterial2D::~StructuredMaterial2D() {}

/* It is important to the material averaging algorithm that primitives be stored in an 
 * ordered list according to their layer.  Lower layers are stored first (have priority).
 * This means that once you have added a primitive to a list, you cannot change its
 * layer!
 */
void StructuredMaterial2D::add_primitive(MaterialPrimitive* prim)
{
	std::list<MaterialPrimitive*>::iterator it, insert_pos = _primitives.end();

	if(_primitives.size() == 0) {
		_primitives.push_back(prim);
	}
	else {
		for(it = _primitives.begin(); it != _primitives.end(); it++) {
			if( prim->get_layer() < (*it)->get_layer() ) {
				insert_pos = it; 
				break;
			}
		}
		_primitives.insert(it, prim);
	}
}

void StructuredMaterial2D::add_primitives(std::list<MaterialPrimitive*> primitives)
{
    std::list<MaterialPrimitive*>::iterator it;
    for(it = primitives.begin(); it != primitives.end(); it++) {
        add_primitive(*it);
    }
}

void StructuredMaterial2D::get_values(ArrayXcd& grid, int k1, int k2, int j1, int j2, double sx, double sy)
{
    int N = k2 - k1;

    for(int j = j1; j < j2; j++) {
        for(int k = k1; k < k2; k++) {
            grid((j-j1)*N+k-k1) = get_value(k+sx, j+sy);
        }
    }
}

// This attempts to compute a reasonable average of the materials in a given Yee cell
// Note that there are a few situations where this average will not quite be what they
// should be.  In particular, if three or more materials intersect a cell, this 
// average will begin to deviate from the "correct" average
std::complex<double> StructuredMaterial2D::get_value(double x, double y)
{
	std::complex<double> val = 0.0;
	std::list<MaterialPrimitive*>::iterator it = _primitives.begin();
	MaterialPrimitive* prim;
	GridCell cell;
	
	double xd = x*_dx, //+ _dx/2.0,
		   yd = y*_dy; //+ _dy/2.0;

	double xmin = xd - _dx/2.0,
		   xmax = xd + _dx/2.0,
		   ymin = yd - _dy/2.0,
		   ymax = yd + _dy/2.0,
		   overlap = 1.0;

    bool contains_p1,
         contains_p2,
         contains_p3,
         contains_p4;

	cell.set_vertices(xmin,xmax,ymin,ymax);
	
	if(_primitives.size() == 0) {
		std::cerr << "Error: StructuredMaterial list is empty." << std::endl;
		return 0.0;
	}

	while(it != _primitives.end()) {
		prim = (*it);
        contains_p1 = prim->contains_point(xmin,ymin);
        contains_p2 = prim->contains_point(xmax,ymin);
        contains_p3 = prim->contains_point(xmax,ymax);
        contains_p4 = prim->contains_point(xmin,ymax);
		if(contains_p1 && contains_p2 && contains_p3 && contains_p4 && cell.get_area_ratio() == 1.0)
		{
				return prim->get_material(xd,yd);
		}
		else if(contains_p1 || contains_p2 || contains_p3 || contains_p4)
		{
			overlap = prim->get_cell_overlap(cell);
			val += overlap * prim->get_material(xd,yd);
		}
		it++;
		if(cell.get_area_ratio() == 0) {
			break;
		}
	} // while ends
	// assume background has index of 1.0
	if(cell.get_area_ratio() > 0) {
		val += cell.get_area_ratio()*1.0;
	}
	return val;
}

std::list<MaterialPrimitive*> StructuredMaterial2D::get_primitives()
{
    return _primitives;
}

////////////////////////////////////////////////////////////////////////////////////
// ConstantMaterial2D
////////////////////////////////////////////////////////////////////////////////////
ConstantMaterial2D::ConstantMaterial2D(std::complex<double> value)
{
    _value = value;
}

std::complex<double> ConstantMaterial2D::get_value(double x, double y)
{
    return _value;
}

void ConstantMaterial2D::get_values(ArrayXcd& grid, int k1, int k2, int j1, int j2, double sx, double sy)
{
    int N = k2 - k1;
    for(int i = j1; i < j2; i++) {
        for(int j = k1; j < k2; j++) {
            grid((i-j1)*N + j-k1) = _value;
        }
    }
}

void ConstantMaterial2D::set_material(std::complex<double> val)
{
    _value = val;
}

std::complex<double> ConstantMaterial2D::get_material()
{
    return _value;
}

////////////////////////////////////////////////////////////////////////////////////
// ConstantMaterial3D
////////////////////////////////////////////////////////////////////////////////////
ConstantMaterial3D::ConstantMaterial3D(std::complex<double> value)
{
    _value = value;
}

std::complex<double> ConstantMaterial3D::get_value(double k, double j, double i)
{
    return _value;
}

void ConstantMaterial3D::get_values(ArrayXcd& grid, int k1, int k2, int j1, int j2,
                                    int i1, int i2, double sx, double sy, double sz)
{
    int N = k2 - k1,
        M = j2 - j1;
    for(int i = i1; i < i2; i++) {
        for(int j = j1; j < j2; j++) {
            for(int k = k1; k < k2; k++) {
                grid((i-i1)*N*M + (j-j1)*N + k-k1) = _value;
            }
        }
    }
}

void ConstantMaterial3D::set_material(std::complex<double> val)
{
    _value = val;
}

std::complex<double> ConstantMaterial3D::get_material()
{
    return _value;
}

////////////////////////////////////////////////////////////////////////////////////
// StructuredMaterial3D
////////////////////////////////////////////////////////////////////////////////////

StructuredMaterial3D::StructuredMaterial3D(double X, double Y, double Z,
                                           double dx, double dy, double dz) :
                                           _X(X), _Y(Y), _Z(Z), 
                                           _dx(dx), _dy(dy), _dz(dz)
{
    _background = 1.0;
    _use_cache = true;
    _cache_active = false;
}

// We allocate memory -- Need to free it!
StructuredMaterial3D::~StructuredMaterial3D()
{
	for(auto it = _layers.begin(); it != _layers.end(); it++) {
        delete (*it);
    }
}

void StructuredMaterial3D::add_primitive(MaterialPrimitive* prim, double z1, double z2)
{
    // Dummy variables
    StructuredMaterial2D* layer;
    double znew[2] = {z1, z2},
           z = 0;
    // iterators of the two lists of primitives and Z coordinates
    auto itl = _layers.begin();
    auto itz = _zs.begin();
    // where to insert the primitives and Z coordinates
    std::list<StructuredMaterial2D*>::iterator itl_ins;
    std::list<double>::iterator itz_ins;

    // Make sure the layer has a thickness
    if(z1 == z2) {
        std::cout << "Warning in Structured3DMaterial: Provided layer has no \
                      thickness. It will be ignored." << std :: endl;
        return;
    }
    else if(z2 < z1) {
        std::cout << "Warning in Structured3DMaterial: Provided layer has negative \
                      thickness. It will be ignored." << std :: endl;
        return;
    }
    // If this is the first addition, things are simple
    if(itz == _zs.end()) {
        _zs.push_back(z1);
        _zs.push_back(z2);
        layer = new StructuredMaterial2D(_X, _Y, _dx, _dy);
        layer->add_primitive(prim);
        _layers.push_back(layer);
        return;
    }

    // now we insert the beginning and end point of the layer one at a time, breaking
    // up or inserting new layers as necessary
    for(int i = 0; i < 2; i++) {
        z = znew[i];
        // Primitive's Z at layer with lower bound of itz_ins
        itz = _zs.begin(); itz_ins = _zs.end();
        // Primitive's layer index
        itl = _layers.begin(); itl_ins = _layers.end();
        // figure out where the point is going to go
        while(itz != _zs.end()) {
            if(z >= *itz) { // primitive's Z coord above current layer's bottom
                itz_ins = itz;
                itl_ins = itl;
            }
            itz++;
            if(itl != _layers.end())
                itl++;
        }

        // Three cases to consider: (1) point below stack (2) point above stack (3)
        // point in stack
        if(itz_ins == _zs.end()) {
            layer = new StructuredMaterial2D(_X, _Y, _dx, _dy);
            _layers.push_front(layer);
            _zs.push_front(z);
        }
        else if(itz_ins == --_zs.end() && z != *itz_ins) {
            layer = new StructuredMaterial2D(_X, _Y, _dx, _dy);
            _layers.push_back(layer);
            _zs.push_back(z);
        }
        else {
            // make sure the point to insert is not already in the stack
            if(z != *itz_ins) {
                layer = new StructuredMaterial2D(_X, _Y, _dx, _dy);
                layer->add_primitives( (*itl_ins)->get_primitives() );
                _layers.insert(itl_ins, layer);
                _zs.insert(++itz_ins, z);
            }
        }
    }

    // Finally, insert the supplied MaterialPrimitive into the desired locations
    itz = _zs.begin();
    itl = _layers.begin();

    // figure out where the point is going to go
    while(itl != _layers.end()) {
        z = (*itz);
        if(z >= z1 && z < z2) {
            (*itl)->add_primitive(prim);
        }
        itz++;
        itl++;
    }
}

void StructuredMaterial3D::set_Nsubcell(int Nsubcell)
{
    _Nsubcell = Nsubcell;
}
// Note that this takes a 1D array!
void StructuredMaterial3D::get_values(ArrayXcd& grid, int k1, int k2,
                                                      int j1, int j2,
                                                      int i1, int i2,
                                                      double sx, double sy, double sz)
{
    int Nx = k2-k1, Ny = j2-j1, Nz = i2-i1;
    size_t size = (size_t)Nx*Ny*Nz;
    std::complex<double> *hgrid;
    hgrid = (std::complex<double> *)malloc(size * sizeof(std::complex<double>));

    std::list<double>::iterator itz = _zs.begin(), itz_next;
    auto itl = _layers.begin();

    int primseqwithin, ZLen, LayerLen;
    double *ZList;
    ZLen = _zs.size();
    LayerLen = _layers.size();
    ZList = (double *)malloc(ZLen * sizeof(double));

    int idx = 0;
    for(auto z : _zs){
        ZList[idx] = z;
        idx++;
    }

    std::list<MaterialPrimitive*> primtmp;
    int PrimLen = 0; // total number of primitives across all layers
    int VerticesLen = 0;  // total number of vertices of all primitives
    int *PrimVertNumber, *PrimVertCoordLoc, *PrimLayerSeq, *PrimLayerValue,
                *PrimNumberPerLayer, *PrimIdxPerLayer, *PrimSeqWithinLayer;
    double *PrimVertXAll, *PrimVertYAll;
    double *PrimVertXmax, *PrimVertXmin, *PrimVertYmax, *PrimVertYmin;
    std::complex<double> *PrimMatValue;

    // number of primitives in each layer
    PrimNumberPerLayer = (int *)malloc(LayerLen * sizeof(int));
    PrimIdxPerLayer = (int *)malloc(LayerLen * sizeof(int));

    itl = _layers.begin();
    for(int i=0; i < _layers.size(); i++){
        primtmp = (*itl)->get_primitives();
        PrimIdxPerLayer[i] = PrimLen;
        PrimLen += primtmp.size();
        PrimNumberPerLayer[i] = primtmp.size();
        auto itp = primtmp.begin();
        for(int j=0; j < primtmp.size(); j++){
            VerticesLen += (*itp)->get_number_vertices();
            itp++;
        }
        itl++;
    }

    // number of vertices for each primitive
    PrimVertNumber = (int *)malloc(PrimLen * sizeof(int));
    // location of vertex number in the coordinate arrays
    PrimVertCoordLoc = (int *)malloc(PrimLen * sizeof(int));
    // layer of the corresponding primitive
    PrimLayerSeq = (int *)malloc(PrimLen * sizeof(int));
    // seq no within layer of the corresponding primitive
    PrimSeqWithinLayer = (int *)malloc(PrimLen * sizeof(int));
    // value of the layer def of the corresponding primitive
    PrimLayerValue = (int *)malloc(PrimLen * sizeof(int));
    // material value of ALL primitives
    PrimMatValue = (std::complex<double> *)malloc(PrimLen * sizeof(std::complex<double>));

    // X and Y coordinates of ALL vertices of ALL primitives
    PrimVertXAll = (double *)malloc(VerticesLen * sizeof(double));
    PrimVertYAll = (double *)malloc(VerticesLen * sizeof(double));
    // xmin, xmax, ymin, ymax of ALL primitives
    PrimVertXmin = (double *)malloc(PrimLen * sizeof(double));
    PrimVertXmax = (double *)malloc(PrimLen * sizeof(double));
    PrimVertYmin = (double *)malloc(PrimLen * sizeof(double));
    PrimVertYmax = (double *)malloc(PrimLen * sizeof(double));

    itl = _layers.begin();
    int primidx = 0;
    PrimVertCoordLoc[primidx] = 0;
    for(int i=0; i < _layers.size(); i++){
        // list of primitives in the current layer
        primtmp = (*itl)->get_primitives();
        // list head
        auto itp = primtmp.begin();
        // prim seq no within layer
        primseqwithin = 0;
        for(int j=0; j < primtmp.size(); j++){
            PrimVertNumber[primidx] = (*itp)->get_number_vertices();
            PrimLayerValue[primidx] = (*itp)->get_layer();
            PrimLayerSeq[primidx] = i;
            PrimSeqWithinLayer[primidx] = primseqwithin;
            if(primidx==0){ PrimVertCoordLoc[primidx] = 0;}
            else{ PrimVertCoordLoc[primidx] = PrimVertCoordLoc[primidx-1] + PrimVertNumber[primidx-1];}
            PrimMatValue[primidx] = (*itp)->get_material(0.0,0.0);

            (*itp)->get_vertices(PrimVertXAll+PrimVertCoordLoc[primidx], PrimVertYAll+PrimVertCoordLoc[primidx]);
            PrimVertXmin[primidx] = 1e9; PrimVertXmax[primidx] = -1e9;
            PrimVertYmin[primidx] = 1e9; PrimVertYmax[primidx] = -1e9;
            for(int k=0; k< PrimVertNumber[primidx]; k++){
                if(PrimVertXAll[PrimVertCoordLoc[primidx]+k]<PrimVertXmin[primidx]){
                    PrimVertXmin[primidx] = PrimVertXAll[PrimVertCoordLoc[primidx]+k];
                }
                if(PrimVertXAll[PrimVertCoordLoc[primidx]+k]>PrimVertXmax[primidx]){
                    PrimVertXmax[primidx] = PrimVertXAll[PrimVertCoordLoc[primidx]+k];
                }
                if(PrimVertYAll[PrimVertCoordLoc[primidx]+k] < PrimVertYmin[primidx]){
                    PrimVertYmin[primidx] = PrimVertYAll[PrimVertCoordLoc[primidx]+k];
                }
                if(PrimVertYAll[PrimVertCoordLoc[primidx]+k]>PrimVertYmax[primidx]){
                    PrimVertYmax[primidx] = PrimVertYAll[PrimVertCoordLoc[primidx]+k];
                }
            }
            itp++;
            primidx++; primseqwithin++;
        }
        itl++;
    }

    kernelpar *_kpar_host, *_kpar_device;
    _kpar_host = (kernelpar *)malloc(sizeof(kernelpar));
    gpuErrchk(cudaMalloc((void **)&_kpar_device, sizeof(kernelpar)));

    _kpar_host->sz=sz; _kpar_host->sy=sy; _kpar_host->sx=sx;
    _kpar_host->dx=_dx; _kpar_host->dy=_dy; _kpar_host->dz=_dz;
    _kpar_host->i1=i1; _kpar_host->i2=i2; _kpar_host->j1=j1; _kpar_host->j2=j2; _kpar_host->k1=k1; _kpar_host->k2=k2;
    _kpar_host->Nx=Nx; _kpar_host->Ny=Ny; _kpar_host->Nz=Nz; _kpar_host->size = size;
    _kpar_host->PrimLen = PrimLen; _kpar_host->LayerLen = LayerLen;
    _kpar_host->background = _background;

    double *dZList, *dPrimVertXAll, *dPrimVertYAll, *dPrimVertXmin, *dPrimVertXmax, *dPrimVertYmin, *dPrimVertYmax;
    int *dPrimNumberPerLayer, *dPrimIdxPerLayer, *dPrimVertNumber, *dPrimVertCoordLoc, *dPrimLayerSeq, *dPrimSeqWithinLayer, *dPrimLayerValue;
    thrust::complex<double> *dPrimMatValue, *dgrid;
    double *dcellarea, *doverlap;
    bool *dsubcellflag;
    int *dinpoly;

    gpuErrchk(cudaMalloc((void **)&dgrid, size * sizeof(thrust::complex<double>)));
    gpuErrchk(cudaMemset(dgrid, 0.0, size * sizeof(thrust::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dZList, ZLen * sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&dPrimNumberPerLayer, LayerLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimIdxPerLayer, LayerLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertNumber, PrimLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertCoordLoc, PrimLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimLayerSeq, PrimLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimSeqWithinLayer, PrimLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimLayerValue, PrimLen * sizeof(int)));
    gpuErrchk(cudaMalloc((void **)&dPrimMatValue, PrimLen * sizeof(thrust::complex<double>)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertXAll, VerticesLen * sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertYAll, VerticesLen * sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertXmin, PrimLen * sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertXmax, PrimLen * sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertYmin, PrimLen * sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&dPrimVertYmax, PrimLen * sizeof(double)));

    gpuErrchk(cudaMemcpy(dZList, ZList, ZLen * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimNumberPerLayer, PrimNumberPerLayer, LayerLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimIdxPerLayer, PrimIdxPerLayer, LayerLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertNumber, PrimVertNumber, PrimLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertCoordLoc, PrimVertCoordLoc, PrimLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimLayerSeq, PrimLayerSeq, PrimLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimSeqWithinLayer, PrimSeqWithinLayer, PrimLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimLayerValue, PrimLayerValue, PrimLen * sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimMatValue, PrimMatValue, PrimLen * sizeof(thrust::complex<double>), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertXAll, PrimVertXAll, VerticesLen * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertYAll, PrimVertYAll, VerticesLen * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertXmin, PrimVertXmin, PrimLen * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertXmax, PrimVertXmax, PrimLen * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertYmin, PrimVertYmin, PrimLen * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dPrimVertYmax, PrimVertYmax, PrimLen * sizeof(double), cudaMemcpyHostToDevice));

    _kpar_host->grid=dgrid;
    _kpar_host->ZList=dZList;
    _kpar_host->PrimNumberPerLayer = dPrimNumberPerLayer;
    _kpar_host->PrimIdxPerLayer = dPrimIdxPerLayer;
    _kpar_host->PrimVertNumber = dPrimVertNumber;
    _kpar_host->PrimVertCoordLoc = dPrimVertCoordLoc;
    _kpar_host->PrimLayerSeq = dPrimLayerSeq;
    _kpar_host->PrimSeqWithinLayer = dPrimSeqWithinLayer;
    _kpar_host->PrimLayerValue = dPrimLayerValue;
    _kpar_host->PrimMatValue = dPrimMatValue;
    _kpar_host->PrimVertXAll = dPrimVertXAll;
    _kpar_host->PrimVertYAll = dPrimVertYAll;
    _kpar_host->PrimVertXmin = dPrimVertXmin;
    _kpar_host->PrimVertXmax = dPrimVertXmax;
    _kpar_host->PrimVertYmin = dPrimVertYmin;
    _kpar_host->PrimVertYmax = dPrimVertYmax;

    int Nsubcell = _Nsubcell; //128;
    // number of grids within each overlap detection along X axis
    if(Nx==1){ _kpar_host->NcellBlk = 1; }
    else{ _kpar_host->NcellBlk = 4; }
    _kpar_host->Nsubcell = Nsubcell;  // number of subcell points on X/Y
    // preallocation of 2D temporary arrays on the XY plane
    gpuErrchk(cudaMalloc((void **)&dsubcellflag, (size_t)Nx*Ny*Nsubcell*Nsubcell*sizeof(bool)));
    gpuErrchk(cudaMalloc((void **)&dcellarea, size*sizeof(double)));
    gpuErrchk(cudaMalloc((void **)&doverlap, size*sizeof(double)));
    _kpar_host->subcellflag = dsubcellflag;
    _kpar_host->cellarea = dcellarea;
    _kpar_host->overlap = doverlap;

    gpuErrchk(cudaMalloc((void **)&dinpoly, Nx*Ny*sizeof(int)));
    _kpar_host->inpoly = dinpoly;

    cudaError_t cudaerr;
    gpuErrchk(cudaMemcpy(_kpar_device, _kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));

    cudaDeviceSetLimit(cudaLimitMallocHeapSize, 256*1024*1024);

    // laminated launch
    dim3 block_dim(8, 4);
    dim3 grid_dim((int)ceil(Nx/8.0 / _kpar_host->NcellBlk), (int)ceil(Ny/4.0), Nz);
    _kpar_host->Nz = 1;
    for(int i=i1; i<i2; i++){
        _kpar_host->i2 = i;
        gpuErrchk(cudaMemcpy(_kpar_device, _kpar_host, sizeof(kernelpar), cudaMemcpyHostToDevice));
//        kernel_get_values_Ncell <<< grid_dim, block_dim >>> (_kpar_device);
        invoke_kernel_get_values_Ncell(_kpar_device, grid_dim, block_dim);
        cudaerr = cudaDeviceSynchronize();
        if(cudaerr!=cudaSuccess){ printf("kernel launch status \"%s\".\n", cudaGetErrorString(cudaerr)); }
    }
//     cudaerr = cudaGetLastError();
//     printf("last status \"%s\".\n", cudaGetErrorString(cudaerr));

    gpuErrchk(cudaMemcpy(grid.data(), dgrid, size*sizeof(thrust::complex<double>), cudaMemcpyDeviceToHost));

    gpuErrchk(cudaFree(dgrid));
    gpuErrchk(cudaFree(dZList));
    gpuErrchk(cudaFree(dPrimNumberPerLayer));
    gpuErrchk(cudaFree(dPrimIdxPerLayer));
    gpuErrchk(cudaFree(dPrimVertNumber));
    gpuErrchk(cudaFree(dPrimVertCoordLoc));
    gpuErrchk(cudaFree(dPrimLayerSeq));
    gpuErrchk(cudaFree(dPrimSeqWithinLayer));
    gpuErrchk(cudaFree(dPrimLayerValue));
    gpuErrchk(cudaFree(dPrimMatValue));
    gpuErrchk(cudaFree(dPrimVertXAll));
    gpuErrchk(cudaFree(dPrimVertYAll));
    gpuErrchk(cudaFree(dPrimVertXmin));
    gpuErrchk(cudaFree(dPrimVertXmax));
    gpuErrchk(cudaFree(dPrimVertYmin));
    gpuErrchk(cudaFree(dPrimVertYmax));
    gpuErrchk(cudaFree(dsubcellflag));
    gpuErrchk(cudaFree(dcellarea));
    gpuErrchk(cudaFree(doverlap));
    gpuErrchk(cudaFree(dinpoly));
    gpuErrchk(cudaFree(_kpar_device));

    delete[] ZList; delete[] PrimNumberPerLayer; delete[] PrimIdxPerLayer;
    delete[] PrimVertNumber; delete[] PrimVertCoordLoc; delete[] PrimLayerSeq; delete[] PrimSeqWithinLayer;
    delete[] PrimLayerValue; delete[] PrimMatValue; delete[] PrimVertXAll; delete[] PrimVertYAll;
    delete[] PrimVertXmin; delete[] PrimVertXmax; delete[] PrimVertYmin; delete[] PrimVertYmax;
    delete[] _kpar_host;
}


std::complex<double> StructuredMaterial3D::get_value(double k, double j, double i)
{
    double zmin = (i-0.5) * _dz, zmax = (i+0.5) * _dz;

    std::complex<double> value = 0.0, mat_val;
    std::list<double>::iterator itz = _zs.begin(),
                                itz_next;
    auto itl = _layers.begin();
    auto itcv = _cached_values.begin();
    auto itcf = _cached_flags.begin();
    bool cached = false;
    int jc = 0, kc = 0;
    // Check if i is below the stack
    if(zmax <= *itz) { return _background; }
    if(zmax > *itz && zmin < *itz) {
        value = (*itz - zmin) / _dz * 1.0;
        zmin = *itz;
    }
    while(itl != _layers.end())
    {
        itz_next = std::next(itz);
        if(zmin >= *itz && zmax <= *itz_next)
        {
            if(_use_cache and _cache_active) {
                jc = int(j) - _cache_j0;
                kc = int(k) - _cache_k0;
                cached = (*itcf)(jc, kc);
                if(cached) {
                    mat_val = (*itcv)(jc, kc);
                }
                else {
                    mat_val = (*itl)->get_value(k, j);
                    (*itcv)(jc, kc) = mat_val;
                    (*itcf)(jc, kc) = true;
                }
            }
            else {
                mat_val = (*itl)->get_value(k, j);
            }
            value += (zmax - zmin) / _dz * mat_val;
            return value;
        }
        else if(zmin >= *itz && zmin < *itz_next && zmax > *itz_next)
        {
            if(_use_cache and _cache_active) {
                jc = int(j) - _cache_j0;
                kc = int(k) - _cache_k0;

                cached = (*itcf)(jc, kc);
                if(cached) {
                    mat_val = (*itcv)(jc, kc);
                }
                else {
                    mat_val = (*itl)->get_value(k, j);
                    (*itcv)(jc, kc) = mat_val;
                    (*itcf)(jc, kc) = true;
                }
            }
            else {
                mat_val = (*itl)->get_value(k, j);
            }
            value += (*itz_next - zmin) / _dz * mat_val;
            zmin = *itz_next;
        }
        itl++; itz++; itcv++; itcf++;
    }
    value += (zmax - zmin) / _dz * 1.0;
    return value;
}

/// following taken from Grid_ctypes.cpp

/////////////////////////////////////////////////////////////////////////////////////
// Material2D
/////////////////////////////////////////////////////////////////////////////////////

void Material2D_get_value(Material2D* mat, complex64* val, double x, double y) {
    std::complex<double> value = mat->get_value(x,y);

    val[0].real = std::real(value);
    val[0].imag = std::imag(value);
}

void Material2D_get_values(Material2D* mat, complex64* arr,
        int k1, int k2, int j1, int j2, double sx, double sy)
{
    std::complex<double> val;
    int Ny = j2-j1,
        Nx = k2-k1;

	ArrayXcd grid(Nx*Ny);
    mat->get_values(grid, k1, k2, j1, j2, sx, sy);

    for(int i = 0; i < Nx*Ny; i++) {
        val = grid(i);
        arr[i].real = std::real(val);
        arr[i].imag = std::imag(val);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Grid Material2D
/////////////////////////////////////////////////////////////////////////////////////

GridMaterial2D* GridMaterial2D_new(int M, int N, complex64* arr) {

	ArrayXXcd grid(N,M);
	complex64 val;

	for(int y = 0; y < N; y++) {
		for(int x = 0; x < M; x++) {
			val = arr[y*M+x];
			grid(y,x) = std::complex<double>(val.real, val.imag);
		}
	}

	return new GridMaterial2D(M, N, grid);
}

void GridMaterial2D_delete(GridMaterial2D* mat) {
	delete mat;
}

void GridMaterial2D_set_grid(GridMaterial2D* mat, int M, int N, complex64* arr)
{
	ArrayXXcd grid(N,M);
	complex64 val;

	for(int y = 0; y < N; y++) {
		for(int x = 0; x < M; x++) {
			val = arr[y*M+x];
			grid(y,x) = std::complex<double>(val.real, val.imag);
		}
	}

	mat->set_grid(M, N, grid);
}

int GridMaterial2D_get_M(GridMaterial2D* mat) {
	return mat->get_M();
}

int GridMaterial2D_get_N(GridMaterial2D* mat) {
	return mat->get_N();
}

/////////////////////////////////////////////////////////////////////////////////////
// Structured Material2D
/////////////////////////////////////////////////////////////////////////////////////

StructuredMaterial2D* StructuredMaterial2D_new(double w, double h, double dx, double dy)
{
	return new StructuredMaterial2D(w,h,dx,dy);
}


void StructuredMaterial2D_delete(StructuredMaterial2D* sm)
{
	delete sm;
}


void StructuredMaterial2D_add_primitive(StructuredMaterial2D* sm,
                                      MaterialPrimitive* prim)
{
	sm->add_primitive(prim);
}

/////////////////////////////////////////////////////////////////////////////////////
// MaterialPrimitives
/////////////////////////////////////////////////////////////////////////////////////

void MaterialPrimitive_set_layer(MaterialPrimitive* prim, int layer)
{
	prim->set_layer(layer);
}

int MaterialPrimitive_get_layer(MaterialPrimitive* prim)
{
	return prim->get_layer();
}

void MaterialPrimitive_get_vertices(MaterialPrimitive* prim, double *vx, double *vy)
{
    prim->get_vertices(vx, vy);
}

int MaterialPrimitive_get_number_vertices(MaterialPrimitive* prim)
{
    return prim->get_number_vertices();
}
bool MaterialPrimitive_contains_point(MaterialPrimitive* prim, double x, double y)
{
	return prim->contains_point(x,y);
}

double MaterialPrimitive_get_material_real(MaterialPrimitive* prim,
                                           double x, double y)
{
	return std::real(prim->get_material(x,y));
}

double MaterialPrimitive_get_material_imag(MaterialPrimitive* prim,
                                           double x, double y)
{
	return std::imag(prim->get_material(x,y));
}

/////////////////////////////////////////////////////////////////////////////////////
// Circle Primitives
/////////////////////////////////////////////////////////////////////////////////////

Circle* Circle_new(double x0, double y0, double r) {
	return new Circle(x0, y0, r);
}

void Circle_delete(Circle* c) {
	delete c;
}

void Circle_set_material(Circle* c, double real, double imag)
{
	c->set_material(std::complex<double>(real, imag));
}

void Circle_set_position(Circle* c, double x0, double y0)
{
	c->set_position(x0, y0);
}

void Circle_set_radius(Circle* c, double r)
{
	c->set_radius(r);
}

double Circle_get_x0(Circle* c)
{
	return c->get_x0();
}

double Circle_get_y0(Circle* c)
{
	return c->get_y0();
}

double Circle_get_r(Circle* c)
{
	return c->get_r();
}


/////////////////////////////////////////////////////////////////////////////////////
//Rectangle Primitives
/////////////////////////////////////////////////////////////////////////////////////

Rectangle* Rectangle_new(double x0, double y0, double xspan, double yspan)
{
	return new Rectangle(x0, y0, xspan, yspan);
}

void Rectangle_delete(Rectangle* r) {
	delete r;
}

void Rectangle_set_material(Rectangle* r, double real, double imag)
{
	r->set_material(std::complex<double>(real, imag));
}

void Rectangle_set_position(Rectangle* r, double x0, double y0)
{
	r->set_position(x0, y0);
}

void Rectangle_set_width(Rectangle* r, double width)
{
	r->set_width(width);
}

void Rectangle_set_height(Rectangle* r, double height)
{
	r->set_height(height);
}

/////////////////////////////////////////////////////////////////////////////////////
// Polygon Primitives
/////////////////////////////////////////////////////////////////////////////////////

Polygon* Polygon_new()
{
	return new Polygon();
}

void Polygon_delete(Polygon* poly)
{
	delete poly;
}

void Polygon_add_point(Polygon* poly, double x, double y)
{
	poly->add_point(x,y);
}

void Polygon_add_points(Polygon* poly, double* x, double* y, int n)
{
	poly->add_points(x,y,n);
}


void Polygon_set_point(Polygon* poly, double x, double y, int index)
{
	poly->set_point(x, y, index);
}

void Polygon_set_points(Polygon* poly, double* x, double* y, int n)
{
	poly->set_points(x,y,n);
}

void Polygon_set_material(Polygon* poly, double real, double imag)
{
	poly->set_material(std::complex<double>(real, imag));
}

/////////////////////////////////////////////////////////////////////////////////////
// ConstantMaterial2D
/////////////////////////////////////////////////////////////////////////////////////
ConstantMaterial2D* ConstantMaterial2D_new(double real, double imag)
{
    return new ConstantMaterial2D(std::complex<double>(real, imag));
}

void ConstantMaterial2D_set_material(ConstantMaterial2D* cm, double real, double imag)
{
    cm->set_material(std::complex<double>(real, imag));
}

double ConstantMaterial2D_get_material_real(ConstantMaterial2D* cm)
{
    return std::real(cm->get_material());
}

double ConstantMaterial2D_get_material_imag(ConstantMaterial2D* cm)
{
    return std::imag(cm->get_material());
}

/////////////////////////////////////////////////////////////////////////////////////
// Material3D
/////////////////////////////////////////////////////////////////////////////////////

void Material3D_get_value(Material3D* mat, complex64* val, double x, double y, double z) {
	std::complex<double> value = mat->get_value(x,y,z);

    val[0].real = std::real(value);
    val[0].imag = std::imag(value);
}

void Material3D_get_values(Material3D* mat, complex64* arr, int k1, int k2,
                                                            int j1, int j2,
                                                            int i1, int i2,
                                                            double sx, double sy,
                                                            double sz)
{
    std::complex<double> val;
    int Ny = j2-j1,
        Nx = k2-k1,
        Nz = i2-i1;

	ArrayXcd grid(Nx*Ny*Nz);
    mat->get_values(grid, k1, k2, j1, j2, i1, i2, sx, sy, sz);

    for(int i = 0; i < Nx*Ny*Nz; i++) {
        val = grid(i);
        arr[i].real = std::real(val);
        arr[i].imag = std::imag(val);
    }
}

void Material3D_set_mu(Material3D* mat, complex64* arr, int k1, int k2,
                                                            int j1, int j2,
                                                            int i1, int i2,
                                                            double sx, double sy,
                                                            double sz)
{
    int Ny = j2-j1, Nx = k2-k1, Nz = i2-i1;
    for(int i = 0; i < Nx*Ny*Nz; i++) { arr[i].real = 1.0; }
}

/////////////////////////////////////////////////////////////////////////////////////
// ConstantMaterial3D
/////////////////////////////////////////////////////////////////////////////////////
ConstantMaterial3D* ConstantMaterial3D_new(double real, double imag)
{
    return new ConstantMaterial3D(std::complex<double>(real, imag));
}

void ConstantMaterial3D_set_material(ConstantMaterial3D* cm, double real, double imag)
{
    cm->set_material(std::complex<double>(real, imag));
}

double ConstantMaterial3D_get_material_real(ConstantMaterial3D* cm)
{
    return std::real(cm->get_material());
}

double ConstantMaterial3D_get_material_imag(ConstantMaterial3D* cm)
{
    return std::imag(cm->get_material());
}

/////////////////////////////////////////////////////////////////////////////////////
// Structured3DMaterial
/////////////////////////////////////////////////////////////////////////////////////
StructuredMaterial3D* StructuredMaterial3D_new(double X, double Y, double Z,
                                               double dx, double dy, double dz)
{
    return new StructuredMaterial3D(X, Y, Z, dx, dy ,dz);
}

void StructuredMaterial3D_delete(StructuredMaterial3D* sm)
{
    delete sm;
}

void StructuredMaterial3D_add_primitive(StructuredMaterial3D* sm,
                                        MaterialPrimitive* prim,
                                        double z1, double z2)
{
  sm->add_primitive(prim, z1, z2);
}

void StructuredMaterial3D_set_Nsubcell(StructuredMaterial3D* sm, int Nsubcell)
{
    sm->set_Nsubcell(Nsubcell);
}
/////////////////////////////////////////////////////////////////////////////////////
// MISC
/////////////////////////////////////////////////////////////////////////////////////
void row_wise_A_update(Material2D* eps, Material2D* mu, int ib, int ie, int M, int N,
                       int x1, int x2, int y1, int y2, complex64* vdiag)
{
    int x = 0,
        y = 0,
        j = 0,
        ig = 0,
        Nc = 3,
        component = 0;

    std::complex<double> value;
    std::complex<double> I(0.0, 1.0);

    #pragma omp parallel for private(ig, component, y, x, j)
    for(int i=ib; i < ie; i++) {
        ig = i/Nc;
        component = i-ig*Nc;
        y = ig/N;
        x = ig - y*N;

        j = i-ib;
        if(component == 0) {
            if(x >= x1 && x < x2 && y >= y1 && y < y2) {
                value = I*eps->get_value(x,y);
                vdiag[j].real = std::real(value);
                vdiag[j].imag = std::imag(value);
            }
        }
        else if(component == 1) {
            if(x >= x1 && x < x2 && y >= y1 && y < y2) {
                value = -I*mu->get_value(double(x),double(y)+0.5);
                vdiag[j].real = std::real(value);
                vdiag[j].imag = std::imag(value);
            }
        }
        else {
            if(x >= x1 && x < x2 && y >= y1 && y < y2) {
                value = -I*mu->get_value(double(x)-0.5,double(y));
                vdiag[j].real = std::real(value);
                vdiag[j].imag = std::imag(value);
            }
        }
    }
}
