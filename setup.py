from setuptools import setup

setup(
    name='bleties',
    version='0.1.3',
    description='IES analysis in ciliate genomes with PacBio long reads',
    author='Brandon Seah',
    author_email='kb.seah@tuebingen.mpg.de',
    url='https://github.com/Swart-lab/bleties',
    packages = ['bleties'],
    scripts=['bin/bleties', 'bin/milcor_plot.py', 'bin/milraa_plot.py']
    )
