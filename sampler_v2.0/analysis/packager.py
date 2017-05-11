#!/usr/bin/python

import os
import shutil
import zipfile


exclude_folder_patterns = ["logs", 'build', 'data']
include_folders = ["sampler_v2.0", "clustering", "data", "distance", "src"]
include_file_types = ["cpp", "sh", "py", "c", "h", "r", "cc"]
include_file_names = ["Makefile"]
include_file_patterns = ["compiled_"]

fileList = []


def filename_matches_list(filename):
    global include_file_names

    if filename in include_file_names:
        return True

    for pattern in include_file_patterns:
        if filename.startswith(pattern):
            return True

    return False

def trav(path, depth):
    
    zipf = zipfile.ZipFile("sampler_v2.0/analysis/pkg.zip", "w")
    for dirname, subdirs, files in os.walk(path):
        # print dirname, subdirs, files
        
        sd = False
        for pattern in exclude_folder_patterns:
            if pattern in dirname:
                sd = True
        if sd:
            continue
            
        zipf.write(dirname)
        
        for f in files:
            # zipf.write(os.path.join(dirname, filename))
            if ("." in f and f[f.rindex(".")+1:] in include_file_types) or filename_matches_list(f):
                print os.path.join(dirname, f)
                zipf.write(os.path.join(dirname, f))

    zipf.close()


def main():
    # cwd = os.getcwd()
    # print cwd
    
    # cwd = os.getcwd()
    # print cwd
    

    # return

    # dest = os.path.abspath("..")
    # src = os.path.abspath(".")
    # shutil.copyfile(os.path.join(src, 'looper.sh'),     os.path.join(dest, 'looper.sh'))
    
    os.chdir('../..')
    trav("sampler_v2.0", 0)

    os.chdir('sampler_v2.0')
    # try:
    #     os.remove("looper.sh")
    # except:
    #     print 'could not delete file'

    print "Done"

    

        

    


if __name__ == "__main__":
    main()


    '''
    try:
        os.mkdir("temp_for_copy")
    except:
        pass
    dest = os.path.abspath("temp_for_copy")
    src = os.path.abspath(".")
        

    for f, path in fileList:
        # print os.path.join(cwd, path)," -> ", f
        # print path," -> ", f

        # print os.path.join(dest, path)," -> ", f

        try:
            os.makedirs(os.path.join(dest, path))
        except:
            pass

        # print "copy src:", os.path.join(src, path, f)
        # print "copy dst:", os.path.join(dest, path, f)
        shutil.copyfile(os.path.join(src, path, f),     os.path.join(dest, path, f))


        # print ""
    '''

    