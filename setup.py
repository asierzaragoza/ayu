from setuptools import setup, find_packages

setup(
    name='ayu',
    version='0.5.0',
    url='https://github.com/asierzaragoza/ayu',
    python_requires='>=3.9',
    packages=find_packages(
        include=['ayu', 'data', 'bin']
        ),
    package_data=find_packages(include='data'),
    include_package_data=True,

)