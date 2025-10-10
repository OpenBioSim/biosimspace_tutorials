import requests, os

links = {
    "01": (
        "inputs_01.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/ESjwO8clK9xGg6q__fCfEk4BjoYclwD3xTBeXhd_Cc6iZQ?download=1",
    ),
    "02": (
        "inputs_02.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/ESKt9qodhm9Hr767PMWH6VAB1UinRc_nfWzFTb6vcEYruQ?download=1",
    ),
}


def download(link):
    localfile, url = links[link]
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
