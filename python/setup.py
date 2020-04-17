import os
from numpy.distutils.core import setup, Extension


os.chdir(os.path.dirname(os.path.abspath(__file__)))

srcfiles = ['fermifab_module.c', 'bitfield.c', 'boson_map.c', 'fermi_map.c', 'generate_rdm.c', 'sparse.c', 'tensor_op.c', 'util.c']
module = Extension('fermifab.kernel',
                   sources=['fermifab/src/' + file for file in srcfiles],
                   include_dirs=['fermifab/include', '/usr/include/x86_64-linux-gnu'],
                   extra_link_args=['-lm', '-lblas', '-llapacke'])

setup(
    name='fermifab',
    description='FermiFab toolbox for fermionic many-particle quantum systems',
    version='1.0.0',
    author='Ismael Medina-Suárez, Christian B. Mendl',
    url='https://github.com/cmendl/fermifab',
    packages=['fermifab'],
    ext_modules=[module])
