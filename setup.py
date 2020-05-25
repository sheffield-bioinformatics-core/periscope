from setuptools import setup
from periscope import __version__, _program

setup(
    name='periscope',
    version=__version__,
    packages=['periscope'],
    scripts=['periscope/scripts/Snakefile',
             'periscope/scripts/search_for_sgRNA.py'
             ],
    url='',
    license='',
    author='Matthew Parker',
    author_email='matthew.parker@sheffield.ac.uk',
    description='periscope searches for and quanifies sgRNAs in SARS-CoV-2',
    entry_points="""
    [console_scripts]
    {program} = periscope.periscope:main
    """.format(program=_program),
    include_package_data=True,
)
