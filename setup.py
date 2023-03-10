from setuptools import setup, find_packages
import os

moduleDirectory = os.path.dirname(os.path.realpath(__file__))


def readme():
    with open(moduleDirectory + '/README.md') as f:
        return f.read()


setup(
    name="gkligo",
    description='GW (LIGO) utilities',
    long_description=readme(),
    long_description_content_type="text/markdown",
    version="0.0.2",
    author='genghisken',
    author_email='ken.w.smith@gmail.com',
    license='MIT',
    url='https://github.com/genghisken/gkligo',
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.9',
          'Topic :: Utilities',
    ],
    install_requires=[
          'gcn-kafka',
          'healpy',
          'gkutils',
          'docopt',
          'numpy',
          'ligo.skymap',
          'astropy',
          'mocpy',
          'pyYAML',
      ],
    python_requires='>=3.7',
    entry_points = {
        'console_scripts': ['downloadGWAlerts=gkligo.scripts.python.downloadGWAlerts:main'],
    },
)
