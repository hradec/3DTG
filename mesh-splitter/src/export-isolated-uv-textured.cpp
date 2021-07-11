
#include "include/export-isolated-uv-textured.hpp"



using namespace std;
using namespace IECore;
using namespace IECoreImage;
using namespace IECoreScene;
using namespace Imath;


long uv2xy( float u, float v, ImagePrimitivePtr img ){
    long w = long(u * img->getDisplayWindow().max.x);
    long h = long(v * img->getDisplayWindow().max.y);
    return h * img->getDisplayWindow().max.x + w ;
};

string mtl2texture(string mtl_file){
    // find the texture filename inside the mtl file of the obj mesh
    ifstream file(mtl_file);
    for( string line; getline( file, line ); ){
        if(line.find("map_Kd") != string::npos){
            boost::erase_all(line, "map_Kd ");
            if( ! exists(line) ){
                line = dirname(mtl_file)+"/"+line;
            }
            return line;
        }
    }
    return "";
}

int main(int argc, char **argv)
{
    for(int _obj=1;_obj<argc;_obj++){

        string obj     = string(argv[_obj]);
        string mtl     = obj.substr(0, obj.size()-4) + ".mtl";
        string texture = mtl2texture(mtl);
        string outText = texture.substr(0, texture.size()-4) + "-" + basename(obj.substr(0, obj.size()-4)) + ".tif";

        // read obj mesh
        cout << "Loading mesh " << obj << "\n";
        ReaderPtr                 meshReader    = Reader::create( obj );
        MeshPrimitivePtr          mesh          = MeshAlgo::triangulate( meshReader->read().get() );

        // we need to convert float s/t primvars to V2f(UV) primvar,
        // so MeshPrimitiveEvaluator() can work
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

        // and finally we create the new uv primvar!
        // uv is interpolated as FaceVarying
        mesh->variables["uv"] = PrimitiveVariable(PrimitiveVariable::Interpolation::FaceVarying, uvs);

        // now we can finally create the MeshPrimitiveEvaluator
        MeshPrimitiveEvaluatorPtr meshEval      = new MeshPrimitiveEvaluator( mesh );

        // read image texture found in .mtl file
        cout << "Loading texture " << texture << "\n";
        ImageReaderPtr reader   = new ImageReader( texture );
        ImagePrimitivePtr map   = runTimeCast<ImagePrimitive>( reader->read() );
        vector<float> R = map->getChannel<float>( "R" )->readable();
        vector<float> G = map->getChannel<float>( "G" )->readable();
        vector<float> B = map->getChannel<float>( "B" )->readable();
        Box2i size = map->getDataWindow();

        // create a new empty image buffer
        cout << "Creating new empty texture" << "\n";
        ImagePrimitivePtr buffer = new ImagePrimitive( map->getDataWindow(), map->getDisplayWindow() );
        FloatVectorDataPtr bufferR = buffer->createChannel<float> ( "R" );
        FloatVectorDataPtr bufferG = buffer->createChannel<float> ( "G" );
        FloatVectorDataPtr bufferB = buffer->createChannel<float> ( "B" );
        FloatVectorDataPtr bufferA = buffer->createChannel<float> ( "A" );

        // now we rasterize the new empty image buffer, checking
        // each x/y coordinates against the uv set of the mesh,
        // and copy over the RGB data from the original texture
        // only for when the x/y coordinate exists in the mesh UV map
        cout << "Copying over texture data to new texture using UV area..." << "\n";
        PrimitiveEvaluator::ResultPtr  res = meshEval->createResult();
        for( long h=0; h<size.max.y; h++){
            float v = float(h)/size.max.y;
            cout << v << "     \r";
            for( long w=0; w<size.max.x; w++){
                long pindex = h * (size.max.x+1) + w ;
                float u = float(w)/size.max.x;
                if( meshEval->pointAtUV( V2f(u,1-v), res.get()) ){
                    bufferR->writable()[ pindex ] = R[pindex];
                    bufferG->writable()[ pindex ] = G[pindex];
                    bufferB->writable()[ pindex ] = B[pindex];
                    bufferA->writable()[ pindex ] = 1;//R->readable()[ pindex];
                }
            }
        }

        // Matrix filter = getGaussian(3, 3, 10.0);
        // applyFilter(buffer, filter, 4);

        // // dilate
        // int dim=1;
        // for( long h=0; h<size.max.y; h++){
        //     float v = float(h)/size.max.y;
        //     cout << v << "     \r";
        //     for( long w=0; w<size.max.x; w++){
        //         long pindex = h * (size.max.x+1) + w ;
        //         float u = float(w)/size.max.x;
        //         if( bufferA->readable()[ pindex ] == 0  ){
        //             for(int y=dim; y>=-dim; y--){
        //                 for(int x=dim; x>=-dim; x--){
        //                     long _h = min(max(h+y, 0), size.max.y);
        //                     long _w = min(max(w+x, 0), size.max.x);
        //                     long dilate_index = _h * (size.max.x+1) + _w ;
        //                     if( bufferA->readable()[ dilate_index ] > 0  ){
        //                         bufferR->writable()[ pindex ] = bufferR->readable()[ dilate_index ];
        //                         bufferG->writable()[ pindex ] = bufferG->readable()[ dilate_index ];
        //                         bufferB->writable()[ pindex ] = bufferB->readable()[ dilate_index ];
        //                         bufferA->writable()[ pindex ] = 1;
        //                         x=100;
        //                         break;
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }


        // finally, save out our new texture!
        cout << "Writing new texture " << outText << "\n";
        ImageWriterPtr writer = new ImageWriter( buffer,  outText);
        writer->write();
        cout << "Done!"  << "\n\n";

    }
    return 0;
}
