import requests, os

links = {
    "01": (
        "input.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/Ed8O5PzOg5dDkiXleq3cM-8BPWgDN0mQKOvqAJwehEYK5g?download=1",
    ),
    "02": (
        "output.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/ET0OlDt1CAVPgf9X2J59rOcBCXcLZ3IHr1PiKriEaFGj8A?download=1",
    ),
    "03": (
        "example_output.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EbFoczNb57lOqzY9Q2_Us3UB7AkojHxgk2ZPPFxbkqjmNw?download=1",
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
