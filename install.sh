#!/bin/bash

AR2GEMS_PATH=/opt
PYTHON_PY_PATH="${AR2GEMS_PATH}/ar2gems/bin/plugins/Geostat/python/"
PYTHON_UI_PATH="${AR2GEMS_PATH}/ar2gems/bin/plugins/Geostat/"
echo "UI path=${PYTHON_UI_PATH}"
echo "PYTHON path=${PYTHON_PY_PATH}"

cp *.ui $PYTHON_UI_PATH
cp *.py $PYTHON_PY_PATH


