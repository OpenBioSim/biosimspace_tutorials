import requests, os

links = {
    "01": (
        "input.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EcnpFE9iWiNNi1sGkV-ycwQBi42Pea5GPXDe9bEAl-paFw?download=1",
    ),
    "02": (
        "o_xylene_benzene_for_analysis.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EUL_pvV_LfBKql8xljd3BSMBOf2inaY26ME1tgcWhJMKhg?download=1",
    ),
    "03": (
        "exercise_4_5.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/ESvbYVWgophJggIj-1ZwxRsBX-iMOMVGP-CEgDPcsQE0Hw?download=1",
    ),
    "04": (
        "example_output.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EY28QX9JEdlMoeBpZbGE_kcBWtK1Fsj1STFYzoe3lOOpJw?download=1",
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
