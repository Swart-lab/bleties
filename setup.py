from setuptools import setup

# parse version number from bleties/__init__.py
with open('bleties/__init__.py') as fh:
    info = {}
    for line in fh.readlines():
        if line.startswith('__version__'):
            exec(line, info)
            break

setup(
    name='bleties',
    version=info['__version__'],
    description='IES analysis in ciliate genomes with PacBio long reads',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Brandon Seah',
    author_email='kb.seah@tuebingen.mpg.de',
    url='https://github.com/Swart-lab/bleties',
    license='MIT',
    packages = ['bleties'],
    scripts=['bin/bleties', 'bin/milcor_plot.py', 'bin/milraa_plot.py']
    )
