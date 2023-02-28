from typing import Union

import numpy as np
import pandas as pd
from IPython.core.display import display, Image, SVG
from indigo import Indigo
from indigo.renderer import IndigoRenderer
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

indigo = Indigo()
renderer = IndigoRenderer(indigo)

indigo.setOption("render-output-format", "png")
# indigo.setOption("render-superatom-mode", "collapse")
indigo.setOption("render-coloring", True)
# indigo.setOption("render-base-color", "1, 1, 1")
indigo.setOption("render-relative-thickness", 1.5)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def colorize(text, substring, color=bcolors.WARNING):
    return text.replace(substring, f'{color}{substring}{bcolors.ENDC}')


def clear_atom_mapping(rxn):
    for m in rxn.GetReactants():
        for atom in m.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                atom.ClearProp('molAtomMapNumber')
    for m in rxn.GetProducts():
        for atom in m.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                atom.ClearProp('molAtomMapNumber')


def draw_reaction_solvents(df):
    rxn_index = np.random.randint(0, df.shape[0])
    rxn = df.loc[df.index[rxn_index]]
    rxn_mol = ReactionFromSmarts(rxn['reaction_smile'], useSmiles=True)
    clear_atom_mapping(rxn_mol)
    do = rdMolDraw2D.MolDrawOptions()
    do.includeAtomTags = False
    do.includeMetadata = False
    print('Solvent:', rxn['solvents'])
    display(Draw.ReactionToImage(rxn_mol,
                                 subImgSize=(300, 300),
                                 drawOptions=do))
    notes = rxn['notes']
    for (name, _) in rxn['solvents']:
        notes = colorize(notes, name)
    print(notes)


def draw_reaction(data: Union[pd.Series, pd.DataFrame], highlight_text: str = None, render_format='png'):
    if isinstance(data, pd.DataFrame):
        data = data.sample().squeeze()

    rxn = indigo.loadReaction(data['reaction_smile'])
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(rxn)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(rxn)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")

    print("Patent:      ", data['patent'])
    print("Reaction_id: ", data.name)

    if highlight_text:
        print(colorize(data['notes'], highlight_text))
    else:
        print(data['notes'])


def draw_reaction_smi(rxn_smiles: str, render_format='png'):
    rxn = indigo.loadReaction(rxn_smiles)
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(rxn)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(rxn)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")


def draw_mol(mol_smi: str, render_format='png'):
    mol = indigo.loadMolecule(mol_smi)
    indigo.setOption("render-output-format", render_format)
    if render_format == 'png':
        display(Image(renderer.renderToBuffer(mol)))
    elif render_format == 'svg':
        display(SVG(renderer.renderToBuffer(mol)))
    else:
        print(f"{render_format}: unknown render format, 'svg' or 'png' ")
