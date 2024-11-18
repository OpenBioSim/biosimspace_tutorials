import requests, os
import BioSimSpace as BSS

links = {
    "05": (
        "tyk2_atm.tar.bz2",
        os.path.join(BSS.tutorialUrl(), "tyk2_atm.tar.bz2"),
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
