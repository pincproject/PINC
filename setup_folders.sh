#!/bin/bash
#
# Folder Setup Script for mn-fysrp-pic
# Sigvald Marholm, 11.06.15
#
# Description:
#   Builds a typical supporting folder structure around the mn-fysrp-pic directory.
#

rm -rf ../local_data/template
mkdir -p ../local_data/template
cp DiP3D/template/* ../local_data/template/
