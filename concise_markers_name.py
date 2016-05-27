#!/usr/lib/python

def conciseName(mapfile,newfile):
    f = open(mapfile)
    f1 = open(newfile, 'w')
    f1.write(f.readline())
    for i in f:
        j = i.split()
        marker_name = j[0].split('-')
        genotypes = '\t'.join(j[1:])
        pos = len(marker_name)
        if pos > 2:
            newName = marker_name[0]+'-'+marker_name[1]+'-'+marker_name[-1]+'(%s)'%(pos-1)
        else:
            newName = j[0]
        newLines = newName+'\t'+genotypes
        f1.write(newLines+'\n')
    f.close()
    f1.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        conciseName(sys.argv[1],sys.argv[2])
    else:
        print 'Usage:\npython concise_markers_name.py mapfile newfile'
