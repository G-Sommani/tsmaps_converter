import healpy as hp
import numpy as np
import os
import logging
import glob

try: 
    fits_dir = os.environ['NU_FITS_DIR']
except KeyError:
    logging.warning(
        "No NU_FITS_DIR variable set. If you do not set this, importing a "
        "ts map of a neutrino alert will raise an error."
    )

try: 
    numpy_dir = os.environ['NU_NUMPY_DIR']
except KeyError:
    logging.warning(
        "No NU_NUMPY_DIR variable set. If you do not set this, saving a "
        "probability map of a neutrino alert will raise an error."
    )

fits_paths = glob.glob(f"{fits_dir}*.fits.gz")
for path in fits_paths:
    filename = path.split("/")[-1].split(".fits")[0]
    run = int(path.split("/run")[-1].split(".")[0])
    evt = int(path.split(".evt")[-1].split(".")[0])
    llh_values, header = hp.read_map(path, h=True)
    header_dict = dict(header)
    time_mjd = header_dict['EVENTMJD']
    runid = header_dict['RUNID']
    eventid = header_dict['EVENTID']
    ra_deg = header_dict["RA"]  
    dec_deg = header_dict["DEC"]
    mask = ~np.isnan(llh_values)
    probs = np.zeros(len(mask))
    lh_values = np.exp(-llh_values)
    probs[mask] = lh_values[mask]/np.sum(lh_values[mask])
    probs = probs.clip(1.e-16, None)
    n_pixels = float(len(probs))
    nside = hp.pixelfunc.npix2nside(n_pixels)
    index = np.where(probs == np.max(probs))
    (colat, ra_rad) = hp.pix2ang(nside, index)
    dec_rad = np.pi / 2. - colat

    csv_path = "data/IceCube_Gold_Bronze_Tracks.csv"
    try: 
        found = False
        with open(csv_path, 'r') as f:
            file = f.read()
            lines = file.split('\n')[1:-1]
            for line in lines:
                elements = line.split(',')
                csv_runid = int(elements[1])
                csv_eventid = int(elements[2])
                if csv_runid == runid and csv_eventid == eventid:
                    weight = float(elements[14])
                    found = True
        if not found:
            weight = 0.5
            print(f"{runid} {eventid} not found.")
    except KeyError:
        weight = 0.5

    numpy_dict = {
        "RUNID" : runid,
        "EVENTID" : eventid,
        "RA_DEG" : ra_deg,
        "DEC_DEG" : dec_deg, 
        "RA_RAD" : ra_rad,
        "DEC_RAD" : dec_rad,
        "EVENTMJD" : time_mjd,
        "SIGNALNESS" : weight
    }

    np.save(f"{numpy_dir}{filename}_probs", probs)
    np.save(f"{numpy_dir}{filename}_dict", numpy_dict)
