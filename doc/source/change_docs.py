import os

for root,dirnames,filenames in os.walk('docu'):
    for f in filenames:
        if f.endswith('.rst'):
            overwrite = False
            file = open(os.path.join(root,f))
            data = []
            for line in file:
                if ".. automethod:: __init__" in line:
                    overwrite = True
                else:
                    data += [line]
            file.close()
            if overwrite:
                outf = open(os.path.join(root,f),'w')
                for line in data:
                    outf.write(line)
                outf.close()
