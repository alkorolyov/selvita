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
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
chmod +x Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh

# win
set MAMBA_NO_BANNER=1
# linux
export MAMBA_NO_BANNER=1

./mamba init

mamba create -n chem python=3.8 -y
mamba activate chem
mamba install -y numpy pandas dask jupyterlab zstandard colorama matplotlib
pip install ord-schema epam.indigo levenshtein fuzzysearch 

# RXNMApper
# pytorch <1.12 for rxnmapper
mamba install pytorch==1.11.0 torchvision==0.12.0 torchaudio==0.11.0 cudatoolkit=11.3 -c pytorch -y
# mamba install -y rust
pip install rxnmapper
```

