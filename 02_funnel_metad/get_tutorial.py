import requests, os

links = {
    "01": (
        "inputs_01.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/Ed-8fipPSEFMt01VyyZQBBoBpviQb2k5pmzFsdBulaXXUQ?download=1",
    ),
    "02": (
        "inputs_02.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/ES8xuCuWe4dNop5fJeH1opEBoKG0qZZqLph5ZwMadKkJoQ?download=1",
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
