from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='makimono',
      version='0.2.6',
      description='Visualizations for RNA-Seq expression data',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GPL v3 License',
        'Programming Language :: Python :: 2.7',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='GOldwasher GO term enrichment visualization',
      url='http://github.com/hpsbastos/GOldwasher',
      author='hpbastos',
      author_email='hpb29@cam.ac.uk',
      license='GPLv3',
      packages=['GOldwasher'],
      install_requires=[
          'networkx',
          'jinja2',
          'numpy',
          'pandas',
          'matplotlib',
          'rpy2',
          'pydot',
          'graphviz-python'
      ],
      scripts=['bin/makisu'],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      include_package_data=True,
      zip_safe=False)
