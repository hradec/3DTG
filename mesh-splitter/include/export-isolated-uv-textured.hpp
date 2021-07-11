

#include "IECoreImage/ImagePrimitive.h"
#include "IECoreImage/ImageReader.h"
#include "IECoreImage/ImageWriter.h"
#include "IECoreScene/MeshPrimitive.h"
#include "IECoreScene/MeshPrimitiveEvaluator.h"
#include "IECoreScene/MeshAlgo.h"
#include "OpenEXR/ImathBox.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <string>



using namespace std;
using namespace IECore;
using namespace IECoreImage;
using namespace IECoreScene;
using namespace Imath;

long uv2xy( float u, float v, ImagePrimitivePtr img );
string mtl2texture(string mtl_file);


#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

string dirname(const string& s) {
    char sep = '/';
    #ifdef _WIN32
       sep = '\\';
    #endif
    size_t i = s.rfind(sep, s.length());
    if (i != string::npos) {
        return(s.substr(0, i));
    }
    return("");
}

string basename(const string& s) {
    char sep = '/';
    #ifdef _WIN32
       sep = '\\';
    #endif
    size_t i = s.rfind(sep, s.length());
    if (i != string::npos) {
        return(s.substr(i+1));
    }
    return("");
}

bool exists(const string &fileName){
    std::ifstream infile(fileName);
    return infile.good();
}


typedef vector<double> Array;
typedef vector<Array> Matrix;
Matrix getGaussian(int height, int width, double sigma)
{
    Matrix kernel(height, Array(width));
    double sum=0.0;
    int i,j;

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] = exp(-(i*i+j*j)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
            sum += kernel[i][j];
        }
    }

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] /= sum;
        }
    }

    return kernel;
}

void applyFilter(ImagePrimitivePtr image, Matrix &filter){
    int height = image->getDataWindow().max.y;
    int width = image->getDataWindow().max.x;
    int filterHeight = filter.size();
    int filterWidth = filter[0].size();
    int newImageHeight = height-filterHeight+1;
    int newImageWidth = width-filterWidth+1;
    int d,i,j,h,w;

    FloatVectorDataPtr R = image->getChannel<float>( "R" );
    FloatVectorDataPtr G = image->getChannel<float>( "G" );
    FloatVectorDataPtr B = image->getChannel<float>( "B" );
    FloatVectorDataPtr A = image->getChannel<float>( "A" );
    Box2i size = image->getDataWindow();

    ImagePrimitivePtr buffer = new ImagePrimitive( image->getDataWindow(), image->getDisplayWindow() );
    FloatVectorDataPtr bufferR = buffer->createChannel<float> ( "R" );
    FloatVectorDataPtr bufferG = buffer->createChannel<float> ( "G" );
    FloatVectorDataPtr bufferB = buffer->createChannel<float> ( "B" );
    FloatVectorDataPtr bufferA = buffer->createChannel<float> ( "A" );

    for (i=0 ; i<newImageHeight ; i++) {
        cout << float(i)/newImageHeight << "     \r";
        for (j=0 ; j<newImageWidth ; j++) {
            for (h=i ; h<i+filterHeight ; h++) {
                for (w=j ; w<j+filterWidth ; w++) {
                    long pindex = i * (width+1) + j ;
                    long windex = h * (width+1) + w ;
                    bufferR->writable()[pindex] += filter[h-i][w-j]*R->readable()[windex];
                    bufferG->writable()[pindex] += filter[h-i][w-j]*G->readable()[windex];
                    bufferB->writable()[pindex] += filter[h-i][w-j]*B->readable()[windex];
                }
            }
        }
    }

    for( long h=0; h<size.max.y; h++){
        float v = float(h)/size.max.y;
        for( long w=0; w<size.max.x; w++){
            long pindex = h * (size.max.x+1) + w ;
            float u = float(w)/size.max.x;
            R->writable()[ pindex ] = bufferR->readable()[ pindex ];
            G->writable()[ pindex ] = bufferG->readable()[ pindex ];
            B->writable()[ pindex ] = bufferB->readable()[ pindex ];
            A->writable()[ pindex ] = 1;
        }
    }
}

void applyFilter(ImagePrimitivePtr image, Matrix &filter, int times){
    for(int i=0 ; i<times ; i++) {
        applyFilter(image, filter);
    }
}
