#!/bin/sh

sphinx-multiversion docs/source docs/_build
mv docs/_build ./public
cp docs/redirect.html ./public/index.html
