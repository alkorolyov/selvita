# Initial setup steps

### get ord database
sudo apt-get install git
sudo apt-get install git-lfs


git clone https://github.com/open-reaction-database/ord-data
pip install ord-schema

### setup PATH variable
$PATH="$PATH:/workdir/.local/bin"
sudo mcedit ~/.profile
export PATH="$PATH:/workdir/.local/bin"

### install packages
pip install "dask[complete]"

### git clone ORD project
# 
git init
git add *
# create .gitignore
git commit
# create repo in github and attach
git remote add origin https://github.com/username/new_repo
git push -u origin main

# install gh https://github.com/cli/cli/blob/trunk/docs/install_linux.md
gh auth login




### in mambaforge environment
mamba create -n chem python=3.8
mamba activate chem
mamba install -y numpy pandas rdkit dask
pip install ord-schema

