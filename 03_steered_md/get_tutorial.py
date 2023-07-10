import requests,os

links = {"01" : ("data.tar.bz2",
                 "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EQFPq-ebCEtPmtzWyaJrLiwBUjWOBL_QkqquuGa9x8KP5g?download=1"),
        }

def download():
    localfile, url = links["01"]
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
