"""
# Pull example data from internet
# Only needs to be run once to create ./datasets

# https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error
# To get the certificate to work, I ran this:
 cd "/Applications/Python 3.6/"
 sudo "./Install Certificates.command"

# (even though I thought this wasn't the python 3.6 sourced)

#More preferable solution looked like:
import urllib.request as urlrq
import certifi
import ssl
resp = urlrq.urlopen('https://example.com/bar/baz.html',
                      context=ssl.create_default_context(cafile=certifi.where()))
#however urlretrieve is a wrapper for urlopen and doesn't appear to accept the cafile argument
"""

import os
import tarfile
import urllib.request

DOWNLOAD_ROOT = "https://raw.githubusercontent.com/ageron/handson-ml2/master/"
HOUSING_PATH = os.path.join("datasets", "housing")
HOUSING_URL = DOWNLOAD_ROOT + "datasets/housing/housing.tgz"


def fetch_housing_data(housing_url=HOUSING_URL, housing_path=HOUSING_PATH):
    os.makedirs(housing_path, exist_ok=True)
    tgz_path = os.path.join(housing_path, "housing.tgz")
    urllib.request.urlretrieve(housing_url, tgz_path)
    housing_tgz = tarfile.open(tgz_path)
    housing_tgz.extractall(path=housing_path)
    housing_tgz.close()
    return

fetch_housing_data()
