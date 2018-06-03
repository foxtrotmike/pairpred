pymodulespath='/s/chopin/b/grad/minhas/lib64/python/'
pythonpath='/s/chopin/b/grad/minhas/lib/python/'
try:
    import PyML
except ImportError:
    import sys    
    sys.path.append(pymodulespath)
    import matplotlib
    matplotlib.use('Agg') #Must be uncommented for running on department machines
    import os
    if "PYTHONPATH" not in os.environ: #required for yard on department machines
        os.environ["PYTHONPATH"]=pythonpath