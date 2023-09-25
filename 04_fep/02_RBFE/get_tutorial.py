import requests, os

links = {
    "01": (
        "inputs.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EfHshdM9FoVBvAXUp0x1zxMBcoGb4nNcfnSkkxfj51ij6g?download=1",
    ),
    "02": (
        "analysis.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EQ8kWX_hGy5PvzhqvnLWpMwBaqxd_Cd2ez5zEjFI0EfT2g?download=1",
    ),
    "03": (
        "example_output.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EbPiCqbEriZCuhP0SMkRLX0BKdvODAkjIxPDtUf1t8CxmA?download=1",
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
