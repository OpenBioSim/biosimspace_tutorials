import requests, os

links = {
    "01": (
        "input.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EeclcaEfAPlAjtMQzXo3QcIBnZOSX_lX9w81a6yKU6p-NQ?download=1",
    ),
    "02": (
        "o_xylene_benzene_for_analysis.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EYucKLAsmghNknd5914lCTQBC5uQL7fn0Fca_LeOfpaXcA?download=1",
    ),
    "03": (
        "exercise_4_5.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EcD3SVH8VHpMowLj9MvPWksBYWVAlvKLEV_W5g1R3c4n_Q?download=1",
    ),
    "04": (
        "example_output.tar.bz2",
        "https://openbiosim-my.sharepoint.com/:u:/g/personal/director_openbiosim_org/EbXnbu36ozpFq1WNraQeeSQB6wQExM4rFkveOZVkh3dNyw?download=1",
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
