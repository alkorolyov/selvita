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
pip install "dask[complete]"
pip install epam.indigo
pip install pyrobuf
```

git clone ORD project
```
git clone https://github.com/alkorolyov/selvita
```

(optional) in mambaforge environment
```
mamba create -n chem python=3.8
mamba activate chem
mamba install -y numpy pandas rdkit dask
pip install ord-schema
```

