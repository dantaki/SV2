import os
import numpy
import sys
import glob
import fnmatch
from Cython.Build import cythonize
#from distutils.core import setup, Extension
from setuptools import setup,Extension
import distutils.command.install_data
def opj(*args):
        path = os.path.join(*args)
        return os.path.normpath(path)
class wx_smart_install_data(distutils.command.install_data.install_data):
        def run(self):
                install_cmd = self.get_finalized_command('install')
                self.install_dir = getattr(install_cmd, 'install_lib')
                return distutils.command.install_data.install_data.run(self)
def find_data_files(srcdir, *wildcards, **kw):
        OMIT=['.c','.pyc','.egg-info','.so']
        def walk_helper(arg, dirname, files):
                if '.svn' in dirname or 'CVS' in dirname:
                        return
                names = []
                lst, wildcards = arg
                for wc in wildcards:
                        wc_name = opj(dirname, wc)
                        for f in files:
                                filename = opj(dirname, f)
                                if not any(filename.endswith(x) for x in OMIT):
                                        if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                                                names.append(filename)
                if names:
                        lst.append( (dirname, names ) )
        file_list = []
        recursive = kw.get('recursive', True)
        if recursive:
                os.path.walk(srcdir, walk_helper, (file_list, wildcards))
        else:
                walk_helper((file_list, wildcards),
                                        srcdir,
                                        [os.path.basename(f) for f in glob.glob(opj(srcdir, '*'))])
        return file_list
files = find_data_files('sv2/', '*.*')
setup(  name='sv2',
        version='1.4.3.3',
        description='SV2: Support Vector Structural Variation Genotyper',
        url='https://github.com/dantaki/SV2',
        author='Danny Antaki',
        author_email='dantaki@ucsd.edu',
        license='MIT',
        packages=['sv2'],
        package_dir={'sv2': 'sv2/'},
        package_data={ 
            'sv2': ['sv2/resources/*','sv2/resources/training_sets/*',]
            },
        ext_modules=cythonize([
            #Extension('sv2.core',['sv2/core.pyx']),]
            Extension('sv2.FeatureExtraction',['sv2/FeatureExtraction.pyx']),
            Extension('sv2.Genotype',['sv2/Genotype.pyx']),
            Extension('sv2.Preprocess',['sv2/Preprocess.pyx']),
            Extension('sv2.Snv',['sv2/Snv.pyx']),
            Extension('sv2.Svm',['sv2/Svm.pyx']),
            ]
            
        ),
	include_dirs=[numpy.get_include()],
        requires=[
                'cython',
                'json',
		'numpy',
                'pandas',
                'pybedtools',
                'pysam',
		'shutil',
                'wget',
	],
	python_requires='~=2.7',
        include_package_data=True,
	install_requires= [
                'cython',
                'numpy',
		'pandas',
		'pybedtools',
		'pysam>=0.9',
		'scikit-learn>=0.19',
                'wget',
	],
	data_files = files,
	cmdclass = {'install_data': wx_smart_install_data},
        scripts = ['sv2/sv2','sv2/sv2train']
         )
