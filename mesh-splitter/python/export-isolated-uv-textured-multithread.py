#!/usr/bin/env ppython
##########################################################################
#
#  Copyright (c) 2007-2009, Image Engine Design Inc. All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#     * Neither the name of Image Engine Design nor the names of any
#       other contributors to this software may be used to endorse or
#       promote products derived from this software without specific prior
#       written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
#  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
#  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
##########################################################################

# TODO: Not working!

import sys, os
from multiprocessing import Pool, Lock
import imath
import IECore
import IECoreScene
import IECoreImage


class pixel:
    def __init__(self, img):
        self.img = img
        self.w = img.dataWindow.max().x+1
        self.h = img.dataWindow.max().y+1
        self.updated = []

    def xy2index(self, x, y):
        return int((y * self.h * self.w) + (x * (self.w)))
    def get(self, x, y):
        rgb         = [0,0,0]
        pixelIndex  = self.xy2index(x,y)
        rgb = [
            self.img['R'][pixelIndex],
            self.img['G'][pixelIndex],
            self.img['B'][pixelIndex],
        ]
        if 'A' in self.img.keys():
            rgb.append(self.img['A'][pixelIndex])
        return rgb

    def put(self, x, y, rgb):
        pixelIndex  = self.xy2index(x,y)
        self.updated.append(pixelIndex)
        self.img['R'][pixelIndex] = rgb[0]
        self.img['G'][pixelIndex] = rgb[1]
        self.img['B'][pixelIndex] = rgb[2]
        if 'A' in self.img.keys() and len(rgb)>3:
            self.img['A'][pixelIndex] = rgb[3]

    def up2date(self, x, y):
        pixelIndex  = self.xy2index(x,y)
        # print pixelIndex, y, self.w, x
        return  pixelIndex in self.updated

    def puta(self, x, y, a):
        pixelIndex  = self.xy2index(x,y)
        self.img['A'][pixelIndex] = a

def mtl2texture(mtl_file):
    ''' get texture name from mtl file. '''
    if not os.path.exists(mtl_file):
        if not [x for x in sys.argv if '--texture=' in x]:
            Exception( '''
                There is no %s mtl file for the given obj.
                You can specify the texture filename using
                "--texture=<filename>" right after the obj filename.
                for ex:
                    %s %s --texture=objTexture.jpg
            ''' % (mtl_file, os.path.basename(sys.argv[0]), sys.argv[1]))

    textureFilename = [ x.split()[-1].strip('\n').strip() for x in open(mtl_file,'r').readlines() if 'map_Kd' in x ]
    if not textureFilename:
        Exception( 'No texture found in %s file.' % mtl_file)

    ret = textureFilename[0]
    if not os.path.exists(ret):
        ret = os.path.join(os.path.dirname(mtl_file), os.path.basename(textureFilename[0]))
    if not os.path.exists(ret):
        Exception('Cant find texture file %s or %s.' % (textureFilename[0], ret))
    return ret


if __name__ == "__main__":
    if len(sys.argv)<2:
        sys.stderr.write('''
            ERROR: you need to pass the obj file as an argument!
        ''');sys.stdout.flush()
        sys.exit(-1)

    for obj in [x for x in sys.argv[1:] if '--' not in x[0:3]]:
        # grab obj file from command line
        obj = os.path.abspath(obj)
        mtl = obj.replace('.obj','.mtl')

        # grab texture filename from mtl file
        textureName = mtl2texture(mtl)
        outImage = '-'.join([
            os.path.splitext(textureName)[0],
            os.path.basename(os.path.splitext(obj)[0])
        ])+'.tif'

        # load image
        image = IECore.Reader.create( textureName ).read()
        pimage = pixel(image)

        # create a new imagePrimitive for the new image
        outScaledDownBy = 1
        dataWindowResized = imath.Box2i(image.dataWindow.min(),image.dataWindow.max()/outScaledDownBy)
        displayWindowResized = imath.Box2i(image.displayWindow.min(),image.displayWindow.max()/outScaledDownBy)
        tmp=IECoreImage.ImagePrimitive(dataWindowResized, displayWindowResized)
        tmp.createFloatChannel('R')
        tmp.createFloatChannel('G')
        tmp.createFloatChannel('B')
        ptmp = pixel( tmp )

        # load obj file
        mesh = IECore.Reader.create(obj).read()
        mesh = IECoreScene.MeshAlgo.triangulate(mesh)

        # copy separated "s" and "t" primvars to unified V2f "uv" primvar, so
        # we can translate UV coordinates to P using MeshPrimitiveEvaluator
        uv = IECore.V2fVectorData()
        for n in range(0, mesh['s'].data.size()):
        	uv.append(imath.V2f(mesh['s'].data[n], mesh['t'].data[n]))
        mesh["uv"] = IECoreScene.PrimitiveVariable( mesh['s'].interpolation, uv )

        # create a mesh evaluator
        meshEval = IECoreScene.MeshPrimitiveEvaluator( mesh )

        # now we rasterize the new texture by copying
        # the RGB data from the original texture only
        # where MeshPrimitiveEvaluator.pointAtUV()
        # returns a value.
        b = [[0.5,0.5,0.5]]*ptmp.xy2index(1,1)
        print len(b)
        lock = Lock()
        def rasterLine(y,buffer):
            v = float(y) / ptmp.h
            sys.stdout.write("\r%.02f       " % v);sys.stdout.flush()
            for x in range(0, ptmp.w):
                u = float(x) / ptmp.w
                uv = imath.V2f(u, 1.0-v)
                buffer[ptmp.xy2index(u,v)] = [0.5,0.5,0.5]
                # check if the current uv position
                # exists in the uv map of the mesh
                r = meshEval.createResult()
                if meshEval.pointAtUV( uv, r ):
                    # uv returns a mesh P, so
                    # we copy the RGB data!
                    # color = pimage.get(u, v)
                    # ptmp.put( u, v, color )
                    buffer[ptmp.xy2index(u,v)] = pimage.get(u, v)


        p = Pool(5)
        p.map(rasterLine,range(0, ptmp.h), b)
        # for y in range(0, ptmp.h):
        #     rasterLine(y)

        for y in range(0, ptmp.h):
            v = float(y) / ptmp.h
            for x in range(0, ptmp.w):
                u = float(x) / ptmp.w
                sys.stdout.write("\r%.02f   %s %s    " % (v,ptmp.xy2index(u,v),b[ptmp.xy2index(u,v)]));sys.stdout.flush()
                ptmp.put( u, v, b[ptmp.xy2index(u,v)] )


        # now we write out the new created texture.
        sys.stdout.write("Writing new texture %s... " % outImage);sys.stdout.flush()
        IECoreImage.ImageWriter( tmp, outImage ).write()
        sys.stdout.write("Done!\n");sys.stdout.flush()
