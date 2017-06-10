import os, sys

from distutils.core import setup, Extension
from distutils import sysconfig

def parallelCCompile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None, extra_postargs=None, depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N=8 # number of parallel compilations
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects
import distutils.ccompiler
distutils.ccompiler.CCompiler.compile=parallelCCompile

os.chdir('/usr/lib/python3.6/site-packages/DFTools/')
os.environ["CC"] = "/usr/bin/gcc"
#os.environ["CC"] = "ccache clang"
#os.environ["CXX"] = "ccache clang++"
os.environ["CXX"] = "/usr/bin/gcc"
#cpp_args = ['-std=c++14','-stdlib=libc++', '-mmacosx-version-min=10.11']
cpp_args = ['-std=c++14','-fopenmp']
link_args=['-lgomp']
#link_args=[]
ext_modules = [
    Extension(
    'wrap',
        ['calc.cpp','sup_calc.cpp','dip_calc.cpp','ham_calc.cpp', 'wrap.cpp'],
        include_dirs=['pybind11/include','eigen','/usr/include','/usr/lib'],
        library_dirs=['/usr/local/lib'],
        runtime_library_dirs=['/usr/local/lib'],
        libraries=['stdc++'],
    language='c++',
    extra_compile_args = cpp_args,
    extra_link_args = link_args,
    ),
]
import time
t1 = time.time()
setup(
    name='calculations',
    version='0.0.1',
    author='Louis Ponet',
    author_email='louis.ponet@iit.it',
    description='Calculations for GeTe',
    ext_modules=ext_modules,
)
print("building took:{}s".format(time.time()-t1))
