from pathlib import Path
import os
import numpy as np
import pandas as pd
import xarray as xr
from ..config.setup import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def append_to_master(year, month):
    master_path = os.path.join(validation_dir, 'validation_master.nc')
    monthly_path = os.path.join(validation_dir, year,f'validation_{year}_{month}.nc' )



    ds_month = xr.open_dataset(monthly_path)


    #print(ds_month.dims)
    if os.path.exists(master_path):
        ds_master = xr.open_dataset(master_path)
        print(ds_master.dims)
        ds_updated_master = xr.concat([ds_master, ds_month], dim="drifter").sortby("time")
    else:
        print('Your master path does not exist and is assuming that theis is the first time you are running OSCAR')
        ds_updated_master = ds_month

    # make a tmp path with os.path, then atomic replace
    base, ext = os.path.splitext(master_path)      # (".../validation_master", ".nc")
    tmp = f"{base}.tmp{ext}"                       # ".../validation_master.tmp.nc"

    ds_updated_master.to_netcdf(tmp)
    os.replace(tmp, master_path)                   # atomic on POSIX/Windows

    ds_month.close()
    print(f"Updated: {master_path}")


import os
from matplotlib.backends.backend_pdf import PdfPages

import os
from matplotlib.backends.backend_pdf import PdfPages

def save_plots_to_pdf(
    figs,
    explanations,                 # list[str] (same length as figs) OR a single str for all
    pdf_path,
    metadata=None,
    caption_y=0.02,               # bottom margin position (figure coords)
    fontsize=9,
    pad_bottom=0.20,              # extra bottom room for caption
    box=True
):
    """
    Writes a multi-page PDF where each page is a Matplotlib Figure with a caption:
    'Figure N: <explanation>' centered at the bottom. Long captions wrap.

    - figs: list of Matplotlib Figure objects
    - explanations: list of strings (len == len(figs)) or one string used for all
    """
    if isinstance(explanations, str):
        explanations = [explanations] * len(figs)
    if len(explanations) != len(figs):
        raise ValueError("explanations must match number of figs (or be a single string).")

    base, ext = os.path.splitext(pdf_path)
    tmp = f"{base}.tmp{ext}"

    with PdfPages(tmp) as pdf:
        if metadata:
            pdf.infodict().update(metadata)

        for i, (fig, expl) in enumerate(zip(figs, explanations), start=1):
            # Make space for caption (don’t shrink if user already set a larger bottom)
            fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, pad_bottom))

            # Add numbered caption
            text = f"Figure {i}: {expl}"
            bbox = dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.9) if box else None
            caption_artist = fig.text(
                0.5, caption_y, text,
                ha="center", va="bottom",
                fontsize=fontsize,
                wrap=True,
                bbox=bbox,
            )

            pdf.savefig(fig)
            # Clean up so original figs remain unmodified
            caption_artist.remove()

    os.replace(tmp, pdf_path)
