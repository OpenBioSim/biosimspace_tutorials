import requests, os

links = {
    "01": (
        "inputs.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EVlmB0LOdZ9FvFihBzMiUKwBBsHX0HxKQXMM3Si0t1_J7g?download=1",
    ),
    "02": (
        "analysis.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EQJayn0PWv9Nm_Ql8I2XEFUBm7Z4kkPAPDtYLDG0hmmiRQ?download=1",
    ),
    "03": (
        "example_output.tar.bz2",
        "https://openbiosim.sharepoint.com/:u:/s/public/EVCa3MQvEM1Kvc2wTxxZQpgB5ujW1oS3h-c0-OCIS4eGaQ?download=1",
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
