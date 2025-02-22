//////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2007-2010, Image Engine Design Inc. All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of Image Engine Design nor the names of any
//       other contributors to this software may be used to endorse or
//       promote products derived from this software without specific prior
//       written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
//  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//////////////////////////////////////////////////////////////////////////


#include "include/OBJReader2.h"

#include "IECoreScene/MeshPrimitive.h"

#include "IECore/CompoundData.h"
#include "IECore/CompoundParameter.h"
#include "IECore/FileNameParameter.h"
#include "IECore/MessageHandler.h"
#include "IECore/NullObject.h"
#include "IECore/NumericParameter.h"
#include "IECore/ObjectParameter.h"
#include "IECore/SimpleTypedData.h"
#include "IECore/TypedParameter.h"
#include "IECore/VectorTypedData.h"
#include "IECoreScene/MeshAlgo.h"

#include "boost/bind.hpp"
#include "boost/version.hpp"

#if BOOST_VERSION >= 103600
#define BOOST_SPIRIT_USE_OLD_NAMESPACE
#include "boost/spirit/include/classic.hpp"
#else
#include "boost/spirit.hpp"
#endif

#include <fstream>
#include <iostream>
#include <iterator>

using namespace std;
using namespace IECore;
using namespace IECoreScene;
using namespace Imath;
using namespace boost;
using namespace boost::spirit;

IE_CORE_DEFINERUNTIMETYPED(OBJReader2);

// syntactic sugar for specifying our grammar
typedef boost::spirit::rule<boost::spirit::phrase_scanner_t> srule;

const Reader::ReaderDescription<OBJReader2> OBJReader2::m_readerDescription("obj");

OBJReader2::OBJReader2( const std::string &fileName )
	: Reader( "Alias Wavefront OBJ 3D data reader", new ObjectParameter("result", "the loaded 3D object", new
	NullObject, MeshPrimitive::staticTypeId()))
{
	m_fileNameParameter->setTypedValue( fileName );
}

bool OBJReader2::canRead( const string &fileName )
{
	// there really are no magic numbers, .obj is a simple ascii text file

	// so: enforce at least that the file has '.obj' extension
	if(fileName.rfind(".obj") != fileName.length() - 4)
		return false;

	// attempt to open the file
	ifstream in(fileName.c_str());
	return in.is_open();
}

ObjectPtr OBJReader2::doOperation(const CompoundObject * operands)
{
	// for now we are going to retrieve vertex, texture, normal coordinates, faces.
	// later (when we have the primitives), we will handle a larger subset of the
	// OBJ format

	vpf = new IntVectorData();
	m_vpf = &vpf->writable();

	vids = new IntVectorData();
	m_vids = &vids->writable();

	vertices = new V3fVectorData();
	m_vertices = &vertices->writable();

    // texture coordinates as uv
    uv = new V2fVectorData();
    uv->setInterpretation( GeometricData::UV );
    m_uv = &uv->writable();

    // texture coordinates as uv indexed
    IECore::V2fVectorDataPtr uvid = new V2fVectorData();
    uvid->setInterpretation( GeometricData::UV );
    m_uvid = &uvid->writable();
    IECore::IntVectorDataPtr uvids = new IntVectorData();
	m_uvids = &uvids->writable();


    calculate_normals = false;
    // build normals
    V3fVectorDataPtr normals = new V3fVectorData();
	m_normals = &normals->writable();

    // build normals indexed
    V3fVectorDataPtr normalsid = new V3fVectorData();
	m_normalsid = &normalsid->writable();
    IECore::IntVectorDataPtr normalsids = new IntVectorData();
	m_normalsids = &normalsids->writable();


	// parse the file
	parseOBJ();

	// create our MeshPrimitive
	MeshPrimitivePtr mesh = new MeshPrimitive( vpf, vids, "linear", vertices );
    if( uvid->readable().size() )
	{
		mesh->variables.insert(PrimitiveVariableMap::value_type("uv", PrimitiveVariable( PrimitiveVariable::FaceVarying, uv)));
        // cout << "uvid:  " << uvid->readable().size() << endl;
        // cout << "uvids: " << uvids->readable().size() << endl;
		// mesh->variables.insert(PrimitiveVariableMap::value_type("uv", PrimitiveVariable( PrimitiveVariable::FaceVarying, uvid, uvids)));

        if( ! normals->readable().size() )
            m_normals->resize( uv->readable().size() );

	}
    if( normals->readable().size() )
    {
        mesh->variables.insert(PrimitiveVariableMap::value_type("N", PrimitiveVariable(  PrimitiveVariable::FaceVarying, normals)));
    }

    // calculate normals from triangles if no normals in file
    if( calculate_normals )
    {
        cout << endl << "calculateNormals..." << endl;
        MeshAlgo::calculateNormals(mesh.get());
    }

	return mesh;
}


// parse a vertex
void OBJReader2::parseVertex(const char * begin, const char * end)
{
	vector<float> vec;
	srule vertex = "v" >> real_p[append(vec)] >> real_p[append(vec)] >> real_p[append(vec)];
	parse(begin, vertex, space_p);

	// build v
	V3f v;
	v[0] = vec[0];
 	v[1] = vec[1];
 	v[2] = vec[2];

    if( m_vertices->empty() ) cout << "Reading vertices..." << endl;

	// add this vertex
	m_vertices->push_back(v);

    // cout << "\r" << "parseVertex: " << m_vertices->size() << "                ";
}

void OBJReader2::parseVertexColor(const char * begin, const char * end)
{
	vector<float> vec;
	srule vertex = "v" >> real_p[append(vec)] >> real_p[append(vec)] >> real_p[append(vec)] >> real_p[append(vec)] >> real_p[append(vec)] >> real_p[append(vec)];
	parse(begin, vertex, space_p);

	// build v
	V3f v;
	v[0] = vec[0];
 	v[1] = vec[1];
 	v[2] = vec[2];
	V3f rgb;
	rgb[0] = vec[3];
 	rgb[1] = vec[4];
 	rgb[2] = vec[5];

    if( m_vertices->empty() ) cout << "Reading vertices with color..." << endl;

	// add this vertex
	m_vertices->push_back(v);

    // cout << "\r" << "parseVertex: " << m_vertices->size() << "                ";
}

// parse a texture coordinate
void OBJReader2::parseTextureCoordinate(const char * begin, const char * end)
{
	vector<float> vec;
	srule vertex = "vt" >> real_p[append(vec)] >> real_p[append(vec)] >> *(real_p[append(vec)]);
	parse(begin, vertex, space_p);

	// build v
	V3f vt;
	vt[0] = vec[0];
 	vt[1] = vec[1];
	vt[2] = vec.size() == 3 ? vec[2] : 0.0f;

    if( m_introducedTextureCoordinates.empty() ) cout << "Reading Texture Coordinates..." << endl;

	// add this texture coordinate
	m_introducedTextureCoordinates.push_back(vt);
    m_uvid->push_back(V2f(vt[0], vt[1]));

    // cout << "\r" << "parseTextureCoordinate: " << m_introducedTextureCoordinates.size() << "                ";
}

// parse a normal
void OBJReader2::parseNormal(const char * begin, const char * end)
{
	vector<float> vec;
	srule vertex = "vn" >> real_p[append(vec)] >> real_p[append(vec)] >> real_p[append(vec)];
	parse(begin, vertex, space_p);

	// build v
	V3f vn;
	vn[0] = vec[0];
	vn[1] = vec[1];
	vn[2] = vec[2];


    if( m_introducedNormals.empty() ) cout << "Reading Normals..." << endl;

	// add this normal
	m_introducedNormals.push_back(vn);

    // cout << "\r" << "parseNormal: " << m_introducedNormals.size() << "                ";
}

// parse face
void OBJReader2::parseFace(const char * begin, const char * end)
{
	vector<int> vec;
	vector<int> tvec;
	vector<int> nvec;

	srule entry = int_p[append(vec)]
		>> ("/" >> (int_p[append(tvec)] | epsilon_p) | epsilon_p)
        >> ("/" >> (int_p[append(nvec)] | epsilon_p) | epsilon_p);

	srule face = "f"  >> entry >> entry >> entry >> *(entry);
	parse(begin, face, space_p);

        // cout << "\r";
        // cout << vec[0] << " " << vec[1] << " " << vec[2] << " : " << vec.size() << " - ";
        // for(int n=0;n<tvec.size();n++){ cout << tvec[n] << " " ; };
        // cout << " | " << tvec.size() << " | " << nvec.size() << "           ";

    if( m_vpf->empty() ) cout << "Reading Faces..." << endl;

	// push back the degree of the face
	m_vpf->push_back(vec.size());

	// merge in the edges.  we index from 0, so shift them down.
	// also, vertices may be indexed negatively, in which case they are relative to
	// the current set of vertices
	for(vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i)
	{
		m_vids->push_back(*i > 0 ? *i - 1 : m_vertices->size() + *i);
	}

	// merge in texture coordinates and normals, if present
	// OBJ format requires an encoding for faces which uses one of the vertex/texture/normal specifications
	// consistently across the entire face.  eg. we can have all v/vt/vn, or all v//vn, or all v, but not
	// v//vn then v/vt/vn ...
	if(!nvec.empty())
	{
		if(nvec.size() != vec.size())
			throw Exception("invalid face specification");

		// copy in these references to normal vectors to the mesh's normal vector
		for(vector<int>::const_iterator i = nvec.begin(); i != nvec.end(); ++i)
		{
			m_normals->push_back(m_introducedNormals[*i > 0 ? *i - 1 : m_introducedNormals.size() + *i]);
		}
	}
	// otherwise, check if we have specified normals in some previous face
	// if so, and no normals were given here (examples, encoders that do this?), pump in
	// default normal.  the default normal defined here is the zero normal, which is by
	// definition orthogonal to every other vector.  this might result in odd lighting.
	// else
	// {
    //     // if( r_vertices.empty() )
    //     //     r_vertices = vertices->readable();
    //     // V3f l1 = r_vertices[vec[1]] - r_vertices[vec[0]];
    //     // V3f l2 = r_vertices[vec[2]] - r_vertices[vec[0]];
    //     // V3f N = l1 * l2;
    //
    //     calculate_normals = true;
	// 	V3f zero(0.0f, 0.0f, 0.0f);
    //     // zero = N.normalize();
	// 	for(unsigned int i = 0; i < vec.size(); ++i)
	// 	{
    //         // int p0 = i;
    //         // int p1 = i+1 >= vec.size() ? (i+1)-vec.size() : i+1;
    //         // int p2 = i+2 >= vec.size() ? (i+2)-vec.size() : i+2;
    //         // V3f N = (
    //         //     (r_vertices[vec[p1]] - r_vertices[vec[p0]]) *
    //         //     (r_vertices[vec[p2]] - r_vertices[vec[p0]])
    //         // ).normalize();
	// 		m_normals->push_back(zero);
	// 	}
	// }

	//
	// merge in texture coordinates, if present
	//
	if(!tvec.empty())
	{
		if(tvec.size() != vec.size())
			throw Exception("invalid face specification");

		for(unsigned int i = 0; i < tvec.size(); ++i)
		{
			int index = tvec[i] > 0 ? tvec[i] - 1 : m_introducedTextureCoordinates.size() + tvec[i];
            m_uv->push_back(Imath::V2f(m_introducedTextureCoordinates[index][0], m_introducedTextureCoordinates[index][1]));
            m_uvids->push_back(index);
		}
	}
	else
	{
		for(unsigned int i = 0; i < vec.size(); ++i)
		{
            m_uv->push_back(Imath::V2f(0.0f, 0.0f));
            m_uvids->push_back(0);
		}
	}
}

void OBJReader2::parseGroup(const char *begin, const char *end)
{
	// set current group
	vector<string> groupNames;
	srule grouping = "g" >> *(lexeme_d[alnum_p >> *(alnum_p)][append(groupNames)]);

	parse(begin, grouping, space_p);

	// from 'http://local.wasp.uwa.edu.au/~pbourke/dataformats/obj/':
	// The default group name is default.
	if(groupNames.empty())
	{
		groupNames.push_back("default");
	}

	// \todo associate mesh objects with group names
}

void OBJReader2::parseOBJ() {

	srule comment = comment_p("#");

	// see
	// http://local.wasp.uwa.edu.au/~pbourke/dataformats/obj/

	// vertices
	srule vertex          = ("v"  >> real_p >> real_p >> real_p) [bind(&OBJReader2::parseVertex, this, _1, _2)];
	srule vertex_texture  = ("vt" >> real_p >> real_p)           [bind(&OBJReader2::parseTextureCoordinate, this, _1, _2)];
	srule vertex_normal   = ("vn" >> real_p >> real_p >> real_p) [bind(&OBJReader2::parseNormal,            this, _1, _2)];
	//srule vertex_parameter_space  = "vp" >> real_p >> real_p >> real_p;

	// srule cs_types = ("bmatrix" | "bezier" | "bspline" | "cardinal" | "taylor");
	// srule vertex_curve_or_surface = "cstype" >> "rat" >> cs_types;
	// srule vertex_degree  = "deg" >> real_p >> real_p;
	// srule vertex_basis_matrix  = "bmat";
	// srule vertex_step_size  = "step" >> int_p >> int_p;
	srule vertex_type = vertex | vertex_texture | vertex_normal;


	// vertices with colors
	srule vertex_with_color = ("v"  >> real_p >> real_p >> real_p >> real_p >> real_p >> real_p) [bind(&OBJReader2::parseVertexColor, this, _1, _2)];
	srule vertex_type_color = vertex_with_color | vertex_texture | vertex_normal;


	// elements
	srule point = "p" >> real_p >> *(real_p);
	srule  line = "l" >> int_p >> int_p >> *(int_p);
	srule  face = (ch_p('f') >> *(anychar_p))[bind(&OBJReader2::parseFace, this, _1, _2)];
	// srule curve = "curv";
	// srule curve_2d = "curv2";
	// srule surface = "surf";
	srule element = point | line | face;

 	// free-form curve / surface statements
	// srule parameter = "parm";
	// srule trim_loop = "trim";
	// srule hole_loop = "hole";
	// srule special_curve = "scrv";
	// srule special_point = "sp";
	// srule end_statement = "end";

	// connectivity
	//srule connect = "con";

	// grouping
	srule group_name = ("g" >> *(anychar_p))[bind(&OBJReader2::parseGroup, this, _1, _2)];
	//	srule smoothing_group = "s";
	//	srule merging_group = "mg";
	srule object_name = "o" >> int_p;
	srule grouping = group_name | object_name;

	// display and render attributes
	// srule bevel_interpretation = "bevel";
	// srule color_interpolation = "c_interp";
	// srule dissolve_interpolation = "d_interp";
	// srule level_of_detail = "lod";
	// srule material_name = "usemtl";
	// srule material_library = "mtllib";
	// srule shadow_casting = "shadow_obj";
	// srule ray_tracing = "trace_obj";
	// srule curve_approximation_technique = "ctech";
	// srule surface_approximation_technique = "stech";

	ifstream in(fileName().c_str());
	string str;
	while(getline(in, str))
		parse(str.c_str(), vertex_type_color | vertex_type | element | grouping | comment, space_p);
}
