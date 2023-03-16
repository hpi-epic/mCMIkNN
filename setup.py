import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mCMIkNN",
    version="0.0.1",
    author="Johannes Huegle, Christopher Hagedorn",
    author_email="firstname.lastname@hpi.de",
    description="A library for a conditional independence tests withn constraint-based causal structure learning for mixed discrete and continuous data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JohannesHuegle/ICSL",
    project_urls={
        "Bug Tracker": "https://github.com/JohannesHuegle/ICSL/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.9",
    install_requires = ['networkx==3.0','numpy==1.24.1','pandas==1.5.3','scipy==1.10.0','manm-cs>=0.1.0']
    )
