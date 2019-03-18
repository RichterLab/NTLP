from glob import glob
import os
for f in glob("*.f90"):
    os.rename(f,f.replace(".f90",".F"))
