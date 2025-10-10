import requests, os

links = {
    "01": (
        "data.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EcngMqH4AqRNjsejO0G51lcBO_tCOl3FGmm7Y07M9qUPGw?download=1",
    ),
    "02": (
        "steering.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EfTxzYT2bG9LgSvmEl10r6cBFexFRgKf1572S0IK74TKlg?download=1",
    ),
}


def download(key):
    localfile, url = links[key]
    # Do not download if tarball already found
    if os.path.isfile(localfile):
        return
    print("Downloading %s from openbiosim.org ..." % localfile)
    req = requests.get(url, stream=True)
    with open(localfile, "wb") as f:
        for chunk in req.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
                f.flush()
    print("Extracting compressed tarball ...")
    os.system("tar -xf %s" % localfile)
    # os.system("rm %s" % localfile)
