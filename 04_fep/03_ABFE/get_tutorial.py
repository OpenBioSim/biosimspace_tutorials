import requests,os

links = {"01" : ("input.tar.bz2",
                 "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EZSZkfg60i5OmU1KqvQCSkwBLDu56K9Q_dAreagL3i2IoQ?download=1"),
         "02" : ("output.tar.bz2",
                 "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/Eakue8dBxrVBjBj1p6p6r04BnKhgJaynZ7oCx7UXCut_oA?download=1"),
         "03" : ("example_output.tar.bz2",
                 "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EW0ii-odLAJIqz2exgZxSbMB_hgkfRtSBHmkeH6cnANJLQ?download=1")
        }

def download(link):
    localfile, url = links[link]
    # Do not download if tarball already found
    if os.path.isfile(localfile):
        return 
    print ("Downloading %s from openbiosim.org ..." % localfile)
    req = requests.get(url, stream=True)
    with open(localfile, 'wb') as f:
        for chunk in req.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
                f.flush()
    print ("Extracting compressed tarball ...")
    os.system("tar -xf %s" % localfile)
    #os.system("rm %s" % localfile)
