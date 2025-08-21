# OSCAR Setup Instructions

## Required Python Modules
Here is a shortlist of the packages to install. 

xarray , pyyaml,scipy, copernicusmarine , numpy , pandas, cdsapi, matplotlib ,cartopy , scikit-learn, podaac-data-subscriber  


---

## Accounts & Authentication

### 1. NASA Earthdata (PO.DAAC)
- Create an account at [Earthdata Login](https://urs.earthdata.nasa.gov/).
- Create a `.netrc` file in your home directory:

```bash
nano ~/.netrc
```

Add the following (replace with your own credentials):

```
machine urs.earthdata.nasa.gov
  login YOUR_USERNAME
  password YOUR_PASSWORD
```

---

### 2. Copernicus Marine
- Register at [Copernicus Marine](https://marine.copernicus.eu/).
- Add your credentials to the same `.netrc` file:

```
machine auth.marine.copernicus.eu
  login YOUR_USERNAME
  password YOUR_PASSWORD
```

---

### 3. Copernicus Climate Data Store (CDS)
- Create an account at the [CDS Climate Data Store](https://cds.climate.copernicus.eu/).
- Accept the dataset licence(s) under the **Licence** tab.
- Copy your **personal access token** from your [profile page](https://cds.climate.copernicus.eu/how-to-api).
- Create a `.cdsapirc` file:

```bash
nano ~/.cdsapirc
```

Add:

```
url: https://cds.climate.copernicus.eu/api
key: <YOUR-PERSONAL-ACCESS-TOKEN>
```

---

## macOS Note
If you are on macOS and encounter SSL certificate issues, run:

```bash
open "/Applications/Python 3.12/Install Certificates.command"
```

---

## Quick Start
Once everything is set up, run the package as a module from the project root:

```bash
python3 -m pyoscar.main
```
