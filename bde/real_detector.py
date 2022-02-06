# Several molecules in the existing compound database are not physically realistic - eg. O_6.
# To detect them, lets use api calls for PubChem - if there's experimental data, go ahead with compound
# processing. We'll use the PubChem View REST API (the normal PUG REST API doesn't allow for compound
# experimental data to be retrieved) and look for a 404 Not Found API error to screen.
#
# Sabari Kumar, 2/6/22
#
#
import requests
import logging

PUBCHEM_VIEW_BASE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'
PUBCHEM_VIEW_BASE_URL_2 = '/JSON?heading=Experimental+Properties'
PUBCHEM_CIDS_BASE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'
PUBCHEM_CIDS_BASE_URL_2 = '/cids/TXT'

logging.getLogger("cclib").setLevel(30)


def getpubchemid(smiles=None):
    # Returns a PubChem numerical CID identifier for a given smiles string. Used to compose the PubChem API request
    # URL. As of Feb 2022, The PubChem API returns a byte-encoded CID string with a newline character at the end.
    # Params:
    #   smiles : (string) Input SMILES string
    # Returns:
    #   pubchemid: (int) PubChem numerical CID.
    try:
        response = requests.get(PUBCHEM_CIDS_BASE_URL + smiles + PUBCHEM_CIDS_BASE_URL_2)
        response.raise_for_status()
        return response.content.strip().decode('utf-8')
    except requests.exceptions.HTTPError as rh:
        logging.error('Requested CID not found in PubChemDB. Request SMILES: ' + smiles)


def checkpubchemid(pubchemid=None):
    try:
        response = requests.get(PUBCHEM_VIEW_BASE_URL + str(pubchemid) + PUBCHEM_VIEW_BASE_URL_2)
        response.raise_for_status()
        return True
    except requests.exceptions.HTTPError:
        return False
    except requests.exceptions.ConnectionError:
        logging.error('Cannot connect to PubChem DB: Skipping validation step. Request CID: ' + str(pubchemid))
