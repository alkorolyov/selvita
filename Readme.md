# Initial setup steps

get ord database
```
sudo apt-get install git
sudo apt-get install git-lfs

git clone https://github.com/open-reaction-database/ord-data
pip install ord-schema
```

setup PATH variable
```
$PATH="$PATH:/workdir/.local/bin"
sudo mcedit ~/.profile
export PATH="$PATH:/workdir/.local/bin"
```

install packages
```
pip install "dask[complete]" epam.indigo zstandard colorama matplotlib
pip install levenshtein fuzzysearch
```

git clone ORD project
```
git clone https://github.com/alkorolyov/selvita
```

(optional) in mambaforge environment
```
set MAMBA_NO_BANNER=1
mamba create -n chem
mamba activate chem
mamba install -y numpy pandas dask jupyterlab zstandard colorama matplotlib
pip install ord-schema epam.indigo
pip install levenshtein fuzzysearch

```

