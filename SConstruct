import os
import shutil

env=Environment(tools = ['mingw'], lang = ['english']);

#env['CUDA_TOOLKIT_PATH'] = '/opt/cuda'
#env.Tool('cuda', toolpath=['tools'])

env.Append(
    CCFLAGS     = '-O2 -fopenmp -fPIC',
    CXXFLAGS    = '--std=gnu++0x -fopenmp -fPIC -g -Wall',	#-shared
    #NVCCFLAGS   = '-arch=sm_20 -m64 -Xcompiler -fopenmp -Xcompiler -fPIC -I../dep/include',
    CPPPATH		= ['E:\opencv\GC32243\include', 'E:\gdal\include'],
    LIBPATH		= ['E:/mingw/bin', 'E:\\opencv\\GC32243\\bin', 'E:\gdal\GC32'],
    LIBS		= ['gdal-1'],
    LINKFLAGS	= ['-fopenmp -fPIC']
)

#env.Library('rdacc', Glob('../src/*.cpp') + Glob('../src/*.cu') + Glob('../src/*/*.cpp') + Glob('../src/*/*.cu') + Glob('../src/*/*/*.cu') + Glob('../src/*/*/*.cpp'))
#shutil.copyfile('librdacc.a', '../bin/librdacc.a')

env.Program('defog', Glob('src/*.cpp'))
