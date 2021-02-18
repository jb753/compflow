#!/bin/bash
# Deploy the site to the SRFC
rsync -r _build/html/ jb753@srcf:~/public_html/compflow-docs/
