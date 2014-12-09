from distutils.core import setup

setup(
    name='bfieldsim',
    version='0.1dev',
    packages=['bfieldsim',],
    license='MIT',
    long_description=open('README.txt').read(),
    author='Richard W. Turner',
    author_email='rwturner@stanford.edu',
    install_requires = ['numdifftools', 'numpy', 'matplotlib', 'scipy'],
    description='Code for simulating fields due to atom chips'
)
