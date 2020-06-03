import os
import sys

# Add the package root directory to PYTHONPATH
# and import a copy of our module so tests can use it
PKG_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PKG_ROOT)
import compflow

TEST_DIR = os.path.join(PKG_ROOT,'test')
