from setuptools import setup, find_packages

setup(
    name='dotate',
    version='1.2.0',
    description='A tool for annotating protein domains based on HMMsearch domain-table output.',
    author='Andre Lecona Buttelli',
    author_email='andrelecona@elsi.com',
    url='https://github.com/AndreLecona/Dotate',
    license="GPL-3.0-or-later",
    packages=find_packages(include=['dotate_core', 'dotate_core.*']),
    install_requires=[
        'pandas>=1.0',
        'tqdm>=4.0',
        'sqlalchemy>=1.3',
    ],
    entry_points={
        'console_scripts': [
            'dotate=dotate_core.argparse:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
    python_requires='>=3.6',
)
