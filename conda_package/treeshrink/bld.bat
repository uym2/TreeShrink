R CMD INSTALL -l Rlib dependencies\BMS_0.3.3.tar.gz
if errorlevel 1 exit 1
"%PYTHON%" setup.py install --single-version-externally-managed --record=record.txt
if errorlevel 1 exit 1
