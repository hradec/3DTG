

#include "IECoreImage/ImagePrimitive.h"
#include "IECoreImage/ImageReader.h"
#include "IECoreImage/ImageWriter.h"
#include "IECoreScene/SceneCache.h"
#include "IECoreScene/MeshPrimitive.h"
#include "IECoreScene/MeshPrimitiveEvaluator.h"
#include "IECoreScene/MeshAlgo.h"
#include "IECoreAlembic/AlembicScene.h"
#include "IECoreAlembic/ObjectWriter.h"
#include "IECore/ObjectWriter.h"
#include "IECore/Writer.h"
#include "IECore/BoxOps.h"
#include "IECore/TriangleAlgo.h"
#include "OpenEXR/ImathBox.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <string>

#include "include/OBJReader2.h"
#include "include/OBJReader2_code.h"

/* assimp include files. These three are usually needed. */
#include <assimp/cimport.h>
#include <assimp/scene.h>
#include <assimp/postprocess.h>


using namespace std;
using namespace IECore;
using namespace IECoreImage;
using namespace IECoreScene;
using namespace IECoreAlembic;
using namespace Imath;

#define MAX_VERTEXES 20000

// #define DEBUG


// this macro exists for code clarity only.
// we use the XYZ coordinates as a triple indexed std::map to
// quickly retrieve the vextex index from it's XYZ coordinates.
#define v3fmap(point)           point.x][point.y][point.z

// a cascade of map.find()'s to check if a position already
// exists in the triple indexed std::map of XYZ indexes.
#define contains(_map,item)     ( \
        _map.find(item.x) != _map.end() && \
        _map[item.x].find(item.y) != _map[item.x].end() && \
        _map[item.x][item.y].find(item.z) != _map[item.x][item.y].end() \
)

// a simple macro to check if a value exists in a std::vector
#define vfind(vector,item)      std::find(vector.begin(), vector.end(), item)!=vector.end()



void writeSceneCache( MeshPrimitivePtr mesh, string filename="xx.scc", float time=1.0 ){
    // a snippet function to write Cortex SceneCache files

    // create a new SceneCache Object
    SceneCachePtr root = new SceneCache( filename, IndexedIO::Write );

    // create a transform node using root
    SceneInterfacePtr transformNode = root->child( "transform", SceneInterface::CreateIfMissing );
    M44dDataPtr transformMatrix = new M44dData( M44d().translate( V3d( 0, 0, 0 ) ) );
    transformNode->writeTransform( transformMatrix.get(), time );

    // now create a mesh node, and write out the primitive
    SceneInterfacePtr meshNode = transformNode->child( "mesh", SceneInterface::CreateIfMissing  );
    meshNode->writeObject( mesh.get(), time );
}


void writeAlembicCache( MeshPrimitivePtr mesh, string filename="xx.abc", float time=1.0 ){
    // a snippet function to write Alembic files from Cortex SceneCache objects

    // create a new AlembicScene Object, which essentially is a SceneCache object
    AlembicScene root( filename, IndexedIO::Write );

    // create a transform node using root
    // SceneInterfacePtr transformNode = root.child( "transform", SceneInterface::CreateIfMissing );
    // M44dDataPtr transformMatrix = new M44dData( M44d().translate( V3d( 0, 0, 0 ) ) );
    // transformNode->writeTransform( transformMatrix.get(), time );

    // now create a mesh node, and write out the primitive
    // SceneInterfacePtr meshNode = transformNode->child( "mesh", SceneInterface::CreateIfMissing  );
    SceneInterfacePtr meshNode = root.child( "mesh", SceneInterface::CreateIfMissing  );
    meshNode->writeObject( mesh.get(), time );
}

void UpgradeUVdata( MeshPrimitivePtr mesh ){
    // OBJ reader in Cortex assigns UV coordinates as separated
    // FloatVectorData s and t primitive variables
    // Since MeshPrimitiveEvaluator now user UVs as V2fVectorData primvars,
    // we need to convert s/t primvars to uv.

    // create an empty V2fVectorData, to fill up values and initialize
    // the new uv primvar
    V2fVectorDataPtr uvs = new V2fVectorData();
    uvs->setInterpretation( GeometricData::UV );
    std::vector<Imath::V2f> &uvValues = uvs->writable();

    // acces s and t primvars
    FloatVectorData *s = mesh->variableData<FloatVectorData>( "s", PrimitiveVariable::Interpolation::FaceVarying );
    FloatVectorData *t = mesh->variableData<FloatVectorData>( "t", PrimitiveVariable::Interpolation::FaceVarying );

    // resize uv V2fVectorData to match s size
    uvValues.resize( s->readable().size() );

    // copy s/t values to uv V2fVectorData as V2f elements
    for(std::size_t i = 0; i < s->readable().size(); ++i) {
        uvValues[i] = Imath::V2f( s->readable()[i], t->readable()[i] );
    }

    // and finally we create the new uv primvar
    // uv is interpolated as FaceVarying, just like s/t
    mesh->variables["uv"] = PrimitiveVariable(PrimitiveVariable::Interpolation::FaceVarying, uvs);
}

void splitMesh( MeshPrimitivePtr mesh, vector<MeshPrimitivePtr> &result, vector<string> &resultNames, string resultName){
    // split a mesh primitive in 2 by half bbox at of the major axis

    // if our current mesh has less than the maximum allowed number
    // of vertexes, just return!
    const std::vector<Imath::V3f> &P  = mesh->variableData< V3fVectorData >( "P", PrimitiveVariable::Vertex )->readable();
    if( P.size() <= MAX_VERTEXES ) return;

    // create access variables for input mesh vertex, face and uv data
    const std::vector<int> &vertexIds = mesh->vertexIds()->readable();
    const std::vector<int> &verticesPerFace = mesh->verticesPerFace()->readable();
    // const vector<Imath::V2f> &mesh_uvs = runTimeCast<V2fVectorData>( mesh->variables["uv"].data.get() )->readable();
    const std::vector<Imath::V2f> &mesh_uvs = mesh->expandedVariableData<V2fVectorData>( "uv", PrimitiveVariable::FaceVarying )->readable();
    const std::vector<Imath::V3f> &N  = mesh->expandedVariableData< V3fVectorData >( "N", PrimitiveVariable::FaceVarying )->readable();
    const vector<int> &mesh_uvs_indices = mesh->variables["uv"].indices->readable();
    // const std::vector<Imath::V2f> &uv = mesh->variableData< V2fVectorData >( "uv", PrimitiveVariable::Interpolation::FaceVarying )->readable();

    // Splits the box into two across the major axis, if we don't specify an axis.
    Imath::Box3f halfBbox[2];
    boxSplit( mesh->bound(), halfBbox[0], halfBbox[1] );

#ifdef DEBUG
    cout << "bbox half low: ";
    cout << halfBbox[0].min << " | ";
    cout << halfBbox[0].max << "\n";
    cout << "bbox half high: ";
    cout << halfBbox[1].min << " | ";
    cout << halfBbox[1].max << "\n";
    cout << "\n";
    cout << "P size:" << P.size() << "\n" ;
    cout << "vertex per face size:" << verticesPerFace.size() << "\n" ;
    cout << "uv's size:" << mesh_uvs.size() << " " << mesh_uvs_indices.size() <<  "\n" ;
    cout << "face vertex id size:" << vertexIds.size() << "\n";
    cout << "N size:" << N.size() << "\n" ;
    cout << "\n";
#endif

    // vector<MeshPrimitivePtr> result;
    for( char half=0; half < 2 ; half++ ){

        // faces_vindex stores the initial vertex ID for each face
        // this way we can loop over faces instead of over vertexIDs of faces
        std::vector<int> faces_vindex;
        long vcount=0;
        for( long f=0; f<verticesPerFace.size(); f++){
            faces_vindex.push_back(vcount);
            vcount += verticesPerFace[f];
        }

        // new set of positions for the new mesh to be created
        vector< V3f > p;
        std::vector<int> P2p; // we use to index old P to new p ID
        P2p.resize(P.size());

        // new vertices per face array and vertexIDs array to
        // create the faces for the new mesh
        std::vector<int> vPerFace;
        std::vector<int> vIds;

        // a new set of uv's to fill up with the new uv's for the new mesh
        V2fVectorDataPtr uvData = new V2fVectorData;
    	uvData->setInterpretation( GeometricData::UV );
    	std::vector<Imath::V2f> &uvs = uvData->writable();

        // a new set of normals to fill up with the new normals for the new mesh
        V3fVectorDataPtr nData = new V3fVectorData;
    	nData->setInterpretation( GeometricData::Normal );
    	std::vector<Imath::V3f> &n = nData->writable();

        // mark true for each vertex when they are inside the bbox,
        // and false when outside
        // also create the new positions array for the ones inside (p)
        long _i=0,_o=0;
        std::vector<int> P_is_inside;
        std::map<float, std::map<float, std::map<float, int>>> p_index;
        for( long v=0; v<P.size(); v++){
            bool inside = boxIntersects( halfBbox[half], P[v] );
            P_is_inside.push_back( inside );
            if( inside ){
                _i++;
                p.push_back(  P[v] );
                p_index[v3fmap(P[v])] = p.size()-1;
                P2p[v] = p.size()-1;
            }else{
                _o++;
            }
        }

#ifdef DEBUG
        cout << "ninside: " << _i << " outside: " << _o << "\n";
#endif

        // loop over the faces, checking if each face vertex is inside or
        // outside the bbox.
        for( long f=0; f<faces_vindex.size(); f++){
            int faceNVertex = verticesPerFace[f];
            int fullInside=0;
            vector<int> inside, outside;
            for( char v=0 ; v<faceNVertex ; v++){
                long oldMeshFaceVertexID = faces_vindex[f]+v ;
                if( P_is_inside[ vertexIds[ oldMeshFaceVertexID ] ] ){
                    fullInside++;
                    inside.push_back( oldMeshFaceVertexID );
                }else{
                    outside.push_back( oldMeshFaceVertexID );
                }
            }

            // if fullInside == 0, the face is outside the bbox!
            // so just move on to the next face!
            if( fullInside == 0 ) continue;

            // if fullInside has less than the amount of vertex in the face,
            // some of then are outside the bbox.
            // this means we have to raytrace against the bbox to find the
            // intersection point(s), and create new face(s).
            if( fullInside < faceNVertex ){
                // face is partially inside the bbox
                vector<int> new_edge;
                vector<V3f> new_edge_N;
                vector<V2f> new_edge_uv;
                bool res;

                // trace triangle edges against bbox, and find 2 new vertexes
                // that intersect the bbox plane.
                // the edge needs to be a vector with direction from the
                // ouside vertex to the inside one for this ray-trace to work!
                // destP is the inside vertex
                // origP is the outside vertex
                for(int out=0 ; out<outside.size() ; out++){
                    for(int in=0 ; in<inside.size() ; in++){
                        V3f origP, destP, origN, destN, original_edge, dir, highHitPoint;
                        V2f origUV, destUV;
                        // lets gather all data into variables for clarity!
                        origP  = P[ vertexIds[outside[out]] ];
                        origN  = N[           outside[out]  ];
                        origUV = mesh_uvs[    outside[out]  ];
                        origN  = origN.normalize();
                        destP  = P[ vertexIds[inside[in]] ];
                        destN  = N[           inside[in]  ];
                        destUV = mesh_uvs[    inside[in]  ];
                        destN  = destN.normalize();

                        // and this is the original edge,
                        // with direction from origP vertex to destP
                        // we use it as direction vector for raytrace.
                        dir = original_edge = destP - origP;

                        // passing 3 V3f parameters to boxIntersects triggers
                        // a ray tracing of the first V3f as origina, second as
                        // direction and the last V3f returns the intersection of the
                        // ray with one of the planes of the bbox.
                        // boxIntersects returns True if the ray hit the bbox!
                        res = boxIntersects( halfBbox[half], origP, dir.normalize(), highHitPoint );
                        if(res){
                            // we hit a bbox plane!!
                            // interpolation weight to use for N and UV
                            V3f edge = highHitPoint - origP;
                            float w = edge.length() / original_edge.length();
                            #define linearstep(x,y,w) ((x*(1.0-w))+(y*w))

                            // add a new vertex, if destP is not already in
                            // p array. If probably is since is the vertex inside
                            // the bbox, so we just grab the p array index to
                            // use for face vertex index (p_index)
                            if( ! contains(p_index, destP) ){
                                p.push_back(destP);
                                p_index[v3fmap(destP)] = p.size()-1;
                            }
                            new_edge.push_back(p_index[v3fmap(destP)]);
                            new_edge_N.push_back( destN );
                            new_edge_uv.push_back( destUV );

                            // add the intersection position found by ray
                            // tracing the edge agains the bbox.
                            // use it's position in the p array to
                            // construct new faces
                            // if this position already exists, just grab
                            // the p_index for it!
                            if( ! contains(p_index, highHitPoint) ){
                                p.push_back(highHitPoint);
                                p_index[v3fmap(highHitPoint)] = p.size()-1;
                            }
                            new_edge.push_back(p_index[v3fmap(highHitPoint)]);
                            V3f interpolatedN = linearstep( destN, origN, w).normalize();
                            new_edge_N.push_back( interpolatedN.dot(destN) < 0 ? -interpolatedN : interpolatedN );
                            new_edge_uv.push_back( linearstep( destUV, origUV, w) );

                        }
                    }
                }
                // now we should have 2 edges:
                //    edge 1 - new_edge[0](inside vertex),new_edge[1](traced position)
                //    edge 2 - new_edge[2](inside vertex),new_edge[3](traced position)
                //
                // we always have at least one triangle to add
                // (one vertex inside + 2 intersections)
                vPerFace.push_back( 3 );
                char _vids[3] = {0,1,3};
                for(const char &v : _vids){
                    vIds.push_back( p_index[ v3fmap(p[ new_edge[v] ]) ] );
                    // new interpolated Normals
                    n.push_back( new_edge_N[v] );
                    // new interpolated uv's
                    uvs.push_back( new_edge_uv[v] );
                }

                if( new_edge[0] != new_edge[2] ) {
                    // edge index 0 and 2 are not equal, so we have
                    // 2 vertexes inside + 2 intersections = 4 (quad)
                    // so we have a second triangle to add!!
                    vPerFace.push_back( 3 );
                    char _vids[3] = {0,3,2};
                    for(const char &v : _vids){
                        vIds.push_back( p_index[ v3fmap(p[ new_edge[v] ]) ] );
                        // new interpolated Normals
                        n.push_back( new_edge_N[v] );
                        // new interpolated uv's
                        uvs.push_back( new_edge_uv[v] );
                    }
                }

            }else{
                // face is full inside the bbox
                cout << mesh_uvs[ 0 ] << "\n\n";
                vPerFace.push_back( verticesPerFace[f] );
                for( char v=0 ; v<verticesPerFace[f] ; v++){
                    long oldMeshFaceVertexID = faces_vindex[f]+v ;
                    cout << oldMeshFaceVertexID <<"\n"<< mesh_uvs[ oldMeshFaceVertexID ] << "\n\n";
                    vIds.push_back( P2p[ vertexIds[ oldMeshFaceVertexID ] ] );
                    uvs.push_back( mesh_uvs[ oldMeshFaceVertexID ] );
                    n.push_back( N[ oldMeshFaceVertexID ] );
                }
            }
        }
#ifdef DEBUG
        cout << vIds.size() << "|" << n.size() << "|" << uvs.size() << "\n";
#endif
        // create new meshPrimitive
        MeshPrimitivePtr tmp_mesh = new MeshPrimitive( new IntVectorData(vPerFace), new IntVectorData(vIds), "linear", new V3fVectorData(p) );
        // generate new normals
        // result->variables["N"] = MeshAlgo::calculateNormals(result.get());
        tmp_mesh->variables["N"] = PrimitiveVariable( PrimitiveVariable::FaceVarying, nData );
        // setup new uvs
        tmp_mesh->variables["uv"] = PrimitiveVariable( PrimitiveVariable::FaceVarying, uvData );

#ifdef DEBUG
        cout << result.size() << "\n";
#endif

        char suffix = 65+half;
        string name = resultName + "_" + suffix;
        cout << name << "\n";

        result.push_back( tmp_mesh );
        resultNames.push_back( name );

        splitMesh( tmp_mesh, result, resultNames, name );
    }
}


void MeshPrimitive2GLTF( MeshPrimitivePtr scene )
{
    //
    // code extracted from https://github.com/assimp/assimp/issues/203
    // since there's no documentation on how to do this without
    // importing a pre-existing file that assimp can read.

    // aiMesh *mesh = new aiMesh();
    // mesh->mNumVertices = 3;
    // mesh->mVertices = new aiVector3D [] {{0,0,0}, {0,1,0}, {1,0,0}};
    // mesh->mNumFaces = 1;
    // mesh->mFaces = new aiFace[1];
    // mesh->mFaces[0].mNumIndices = 3;
    // mesh->mFaces[0].mIndices = new unsigned[] { 0, 1, 2 };
    // mesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE; // workaround, issue #3778
    //
    // std::unique_ptr<aiScene> scene(new aiScene());
    // scene->mNumMeshes = 1;
    // scene->mMeshes = new aiMesh * [] { mesh };
    // scene->mNumMaterials = 1;
    // scene->mMaterials = new aiMaterial * [] { new aiMaterial() };
    // scene->mRootNode = new aiNode();
    // scene->mRootNode->mNumMeshes = 1;
    // scene->mRootNode->mMeshes = new unsigned [] { 0 };
    // scene->mMetaData = new aiMetadata(); // workaround, issue #3781

    // create vertices and faces, then pack into an aiMesh
    aiVector3D *vertices = new aiVector3D [] {          // deleted: mesh.h:758
        { -1, -1, 0 },
        { 0, 1, 0 },
        { 1, -1, 0 }
    };

    aiFace *faces = new aiFace[1];                      // deleted: mesh.h:784
    faces[0].mNumIndices = 3;
    faces[0].mIndices = new unsigned [] { 0, 1, 2 };    // deleted: mesh.h:149

    aiMesh *mesh = new aiMesh();                        // deleted: Version.cpp:150
    mesh->mNumVertices = 3;
    mesh->mVertices = vertices;
    mesh->mNumFaces = 1;
    mesh->mFaces = faces;
    mesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE; // workaround, issue #3778

    // ????
    mesh->mTextureCoords[ 0 ] = new aiVector3D[ 3 ];
    mesh->mNumUVComponents[ 0 ] = 3;
    // ????

    // a valid material is needed, even if its empty
    aiMaterial *material = new aiMaterial();            // deleted: Version.cpp:155
    // ???? How to add texture? 

    // a root node with the mesh list is needed; if you have multiple meshes, this must match.
    aiNode *root = new aiNode();                        // deleted: Version.cpp:143
    root->mNumMeshes = 1;
    root->mMeshes = new unsigned [] { 0 };              // deleted: scene.cpp:77

    // pack mesh(es), material, and root node into a new minimal aiScene
    aiScene *out = new aiScene();                       // deleted: by us after use
    out->mNumMeshes = 1;
    out->mMeshes = new aiMesh * [] { mesh };            // deleted: Version.cpp:151
    out->mNumMaterials = 1;
    out->mMaterials = new aiMaterial * [] { material }; // deleted: Version.cpp:158
    out->mRootNode = root;
    out->mMetaData = new aiMetadata(); // workaround, issue #3781

    // and we're good to go. do whatever:
    Assimp::Exporter exporter;
    if (exporter.Export(out, "gltf2", "triangle.gltf") != AI_SUCCESS)
        cerr << exporter.GetErrorString() << endl;
    if (exporter.Export(out, "glb2", "triangle.glb") != AI_SUCCESS)
        cerr << exporter.GetErrorString() << endl;

    // deleting the scene will also take care of the vertices, faces, meshes, materials, nodes, etc.
    delete out;

}



// reference: https://www.youtube.com/watch?v=t6VvtW8y9q4&t=211s&ab_channel=RanjeethMahankali
int main(int argc, char **argv)
{
    for(int _obj=1;_obj<argc;_obj++){

        string obj     = string(argv[_obj]);

        // read obj mesh
        cout << "Loading mesh " << obj << "\n\n";
        // ReaderPtr  meshReader = Reader::create( obj );
        OBJReader2  meshReader( obj );

        // we make sure the geometry is all triangulated here, right after
        // reading the OBJ from disk
        MeshPrimitivePtr mesh = MeshAlgo::triangulate( meshReader.read().get() );

        // Cortex OBJ loader is old, and imports UV data as 2 FloatVectorData primvars named s and t.
        // on latest cortex, all code expects UV data to be stored in a V2fVectorData primvar named uv.
        // this functions "upgrades" the OBJ meshprimtive to the new standard!
        // UpgradeUVdata( mesh );

        // now we can finally create the MeshPrimitiveEvaluator
        MeshPrimitiveEvaluatorPtr meshEval = new MeshPrimitiveEvaluator( mesh );

        // lets start spliting meshes now, and store in a vector.
        vector<string> new_meshes_name;
        vector<MeshPrimitivePtr> new_meshes;
        splitMesh( mesh, new_meshes, new_meshes_name, "chunk" );

#ifdef DEBUG
        cout << new_meshes.size() << "\n";
#endif

        AlembicScene root( "all.abc", IndexedIO::Write );
        for( int n=0; n<new_meshes.size() ; n++){
            // finally, save out our new meshes!
            cout << "Writing new objs " << new_meshes_name[n] << "\n";
            writeAlembicCache( new_meshes[n], new_meshes_name[n]+".abc");
            SceneInterfacePtr meshNode = root.child( new_meshes_name[n], SceneInterface::CreateIfMissing  );

            // save mesh into "all.abc"
            meshNode->writeObject( new_meshes[n].get(), 0 );
        }

        cout << "Done!"  << "\n\n";

    }
    return 0;
}
