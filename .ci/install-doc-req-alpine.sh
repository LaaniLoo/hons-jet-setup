#!/bin/sh

apk update && apk add doxygen git
pip install -U breathe sphinx 
pip install git+https://github.com/Holzhaus/sphinx-multiversion.git
