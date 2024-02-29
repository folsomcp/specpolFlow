from setuptools import setup

setup(
    name='SpecpolFlow',
    version='0.1.0',    
    description='Spectropolarimetric tools',
    url='https://github.com/folsomcp/specpolFlow',
    author='SpecpolFlow Team',
    license='GPL-2.0 license',
    packages=['specpolFlow'],
    install_requires=['numpy',                     
                      'matplotlib',],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',  
    ],
)
