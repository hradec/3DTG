//////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2007-2011, Image Engine Design Inc. All rights reserved.
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

#ifndef IECORESCENE_OBJREADER_H
#define IECORESCENE_OBJREADER_H

#include "IECoreScene/Export.h"
#include "IECoreScene/TypeIds.h"

#include "IECore/Export.h"
#include "IECore/Reader.h"

IECORE_PUSH_DEFAULT_VISIBILITY
#include "OpenEXR/ImathVec.h"
IECORE_POP_DEFAULT_VISIBILITY

#include <vector>

namespace IECoreScene
{

IE_CORE_FORWARDDECLARE(MeshPrimitive);

/// The OBJReader2 class defines a class for reading OBJ mesh data.
/// This is a subset of the full setup of objects encodable in OBJ.
/// \ingroup ioGroup
class IECORESCENE_API OBJReader2 : public IECore::Reader
{
	public:
		IE_CORE_DECLARERUNTIMETYPEDEXTENSION( OBJReader2, 999999, IECore::Reader );

		OBJReader2( const std::string &fileName );

		static bool canRead( const std::string &filename );

	protected:

		IECore::ObjectPtr doOperation( const IECore::CompoundObject * operands) override;

	private:

		static const ReaderDescription<OBJReader2> m_readerDescription;

		// top-level parse method
		void parseOBJ();

		// statement-level parse methods
		void parseVertex(const char *, const char *);
		void parseVertexColor(const char *, const char *);
		void parseTextureCoordinate(const char *, const char *);
		void parseNormal(const char *, const char *);

		// elements
		void parseFace(const char *, const char *);

		// grouping
		void parseGroup(const char *, const char *);

		// pointers to mesh data
        IECore::IntVectorDataPtr vpf;
        IECore::IntVectorDataPtr vids;
        IECore::V3fVectorDataPtr vertices;
        IECore::V2fVectorDataPtr uv;
		std::vector<int> *m_vpf, *m_vids, *m_uvids, *m_normalsids;
		std::vector<Imath::V3f> *m_vertices, *m_normals, *m_normalsid;
		std::vector<float> *m_sTextureCoordinates, *m_tTextureCoordinates;
        std::vector<Imath::V2f> *m_uv, *m_uvid;
        std::vector<Imath::V3f> r_vertices;

        // flag to calculate normals
        bool calculate_normals;

		// local data for assembling mesh
		std::vector<Imath::V3f> m_introducedNormals, m_introducedTextureCoordinates;
};

IE_CORE_DECLAREPTR(OBJReader2);






} // namespace IECoreScene





#endif // IECORESCENE_OBJREADER_H
